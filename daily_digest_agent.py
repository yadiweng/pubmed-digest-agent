import os
import smtplib
from datetime import datetime, timedelta
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

from Bio import Entrez


# =========================================================
# Config (from environment variables / GitHub Secrets)
# =========================================================
EMAIL_SENDER = os.environ.get("EMAIL_SENDER")
EMAIL_PASSWORD = os.environ.get("EMAIL_PASSWORD")
EMAIL_RECEIVER_STR = os.environ.get("EMAIL_RECEIVER")
ENTREZ_EMAIL = os.environ.get("ENTREZ_EMAIL")

DIGEST_TITLE = "Yadi's Daily Biostatistics & Genomics Digest"

# Debug flag: if set to "1", include long queries in email (default OFF)
SHOW_QUERY = os.environ.get("SHOW_QUERY", "0") == "1"

# Optional journal whitelist (comma-separated). If provided, prefer these journals.
# Example value:
#   Nature Methods,Genome Biology,Nature Biotechnology,Bioinformatics,Biostatistics,PNAS,Nature Genetics
JOURNAL_WHITELIST_STR = os.environ.get("JOURNAL_WHITELIST", "")
JOURNAL_WHITELIST = {j.strip().lower() for j in JOURNAL_WHITELIST_STR.split(",") if j.strip()}

# Convert receiver string to list
if EMAIL_RECEIVER_STR:
    EMAIL_RECEIVER = [e.strip() for e in EMAIL_RECEIVER_STR.split(",") if e.strip()]
else:
    EMAIL_RECEIVER = []

# =========================================================
# Topics: DOMAIN ANCHOR TERMS only (to prevent drift)
# =========================================================
TOPICS = {
    "Single-cell statistical methods": [
        "scRNA-seq",
        "single-cell RNA-seq",
        "single cell RNA-seq",
        "single-cell transcriptomics",
        "single cell transcriptomics",
        "single-cell sequencing",
        "single cell sequencing",
    ],
    "Genomics statistical methods": [
        "genome-wide association",
        "GWAS",
        "eQTL",
        "expression quantitative trait locus",
        "fine-mapping",
        "polygenic risk score",
        "polygenic risk",
        "statistical genetics",
    ],
}

# Method-oriented terms (AND condition): general statistical language
METHOD_TERMS = [
    "statistical", "method", "methodology", "model", "modeling",
    "inference", "estimation", "likelihood",
    "bayesian", "posterior", "prior",
    "regression", "penalized", "lasso", "sparse",
    "variational", "probabilistic", "optimization",
    "benchmark", "simulation", "algorithm"
]

# "Method paper" terms (AND condition): signals this is a methods/tool/framework paper
# Include verbs WITHOUT "we" to make matching more robust
METHOD_PAPER_TERMS = [
    "method", "methods", "methodology",
    "propose", "develop", "introduce", "present", "approach",
    "framework", "algorithm", "pipeline", "workflow",
    "software", "tool", "package", "R package", "python package",
    "benchmark", "simulation study", "evaluation", "open-source"
]

# Exclusion terms (NOT condition): reduce disease/clinical discovery papers
# Keep it moderate—too aggressive will cause "no paper found" days.
EXCLUDE_TERMS = [
    "case report", "review", "systematic review", "meta-analysis",
    "patient", "patients", "clinical trial", "randomized",
    "cohort", "prognosis", "survival",
    "glioma", "schizophrenia"
]

CANDIDATE_RETMAX = 200


def _yesterday_pub_date_str() -> str:
    return (datetime.now() - timedelta(days=1)).strftime("%Y/%m/%d")


def _or_clause_titleab_or_mesh(terms: list[str]) -> str:
    clauses = [f'("{t}"[Title/Abstract] OR "{t}"[MeSH Terms])' for t in terms]
    return "(" + " OR ".join(clauses) + ")"


def _or_clause_titleab(terms: list[str]) -> str:
    clauses = [f'"{t}"[Title/Abstract]' for t in terms]
    return "(" + " OR ".join(clauses) + ")"


def _build_query(anchor_terms: list[str], method_terms: list[str], method_paper_terms: list[str],
                 exclude_terms: list[str], pub_date: str) -> str:
    """
    Query structure:
      (anchor terms in Title/Abstract OR MeSH)
      AND (method terms in Title/Abstract)
      AND (method-paper terms in Title/Abstract)
      AND NOT (exclude terms in Title/Abstract)
      AND pub_date[Date - Publication]
    """
    anchor_clause = _or_clause_titleab_or_mesh(anchor_terms)
    method_clause = _or_clause_titleab(method_terms)
    method_paper_clause = _or_clause_titleab(method_paper_terms)

    exclude_clause = ""
    if exclude_terms:
        exclude_clause = " NOT " + _or_clause_titleab(exclude_terms)

    return f"{anchor_clause} AND {method_clause} AND {method_paper_clause}{exclude_clause} AND {pub_date}[Date - Publication]"


def _search_pmids(query: str) -> list[str]:
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=CANDIDATE_RETMAX,
        sort="pub+date"
    )
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])


def _fetch_article_details(pmids: list[str]) -> list[dict]:
    if not pmids:
        return []

    handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="medline", retmode="xml")
    data = Entrez.read(handle)
    handle.close()

    out = []
    for article in data.get("PubmedArticle", []):
        try:
            med = article["MedlineCitation"]["Article"]
            title = med.get("ArticleTitle", "No Title")
            journal = med.get("Journal", {}).get("Title", "No Journal")

            pmid = str(article["MedlineCitation"]["PMID"])
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            author_list = med.get("AuthorList", [])
            if not author_list:
                authors = "No Authors Listed"
            elif len(author_list) == 1:
                a0 = author_list[0]
                authors = f"{a0.get('LastName', '')} {a0.get('Initials', '')}".strip()
            else:
                first = author_list[0]
                last = author_list[-1]
                first_str = f"{first.get('LastName', '')} {first.get('Initials', '')}".strip()
                last_str = f"{last.get('LastName', '')} {last.get('Initials', '')}".strip()
                authors = f"{first_str} ... {last_str}"

            out.append({
                "title": title,
                "authors": authors,
                "journal": journal,
                "journal_lc": journal.lower(),
                "pmid": pmid,
                "url": url,
            })
        except Exception as e:
            print(f"[WARN] Skipping article due to parsing error: {e}")
            continue

    return out


def _select_best_article(articles: list[dict]) -> tuple[dict | None, str]:
    """
    Selection logic:
      1) If JOURNAL_WHITELIST non-empty: pick first article whose journal is in whitelist.
      2) Else: pick first article.
    Returns (article_or_none, note_string).
    """
    if not articles:
        return None, "No paper found."

    if JOURNAL_WHITELIST:
        for a in articles:
            if a["journal_lc"] in JOURNAL_WHITELIST:
                return a, "Preferred journal."
        return articles[0], "Non-preferred journal (fallback)."

    return articles[0], "No journal filter."


def fetch_one_per_topic() -> dict:
    if not ENTREZ_EMAIL:
        print("Error: ENTREZ_EMAIL not set.")
        return {}

    Entrez.email = ENTREZ_EMAIL
    yesterday = _yesterday_pub_date_str()

    results = {}
    for topic, anchor_terms in TOPICS.items():
        query = _build_query(
            anchor_terms=anchor_terms,
            method_terms=METHOD_TERMS,
            method_paper_terms=METHOD_PAPER_TERMS,
            exclude_terms=EXCLUDE_TERMS,
            pub_date=yesterday
        )

        print(f"\n=== Topic: {topic} ===")
        print(f"Query: {query}")

        try:
            pmids = _search_pmids(query)
            articles = _fetch_article_details(pmids)
            picked, note = _select_best_article(articles)
            results[topic] = {"article": picked, "note": note, "query": query}
        except Exception as e:
            print(f"[ERROR] Topic '{topic}' failed: {e}")
            results[topic] = {"article": None, "note": f"Error: {e}", "query": query}

    return results


def format_html_email(topic_payload: dict) -> str:
    yesterday = _yesterday_pub_date_str()

    html = f"""
    <html>
    <head>
    <style>
        body {{ font-family: Arial, sans-serif; }}
        .meta {{ color: #555; margin-bottom: 12px; }}
        table {{ width: 100%; border-collapse: collapse; }}
        th {{ background-color: #1f3b57; color: white; padding: 10px; text-align: left; }}
        td {{ border-bottom: 1px solid #ddd; padding: 8px; vertical-align: top; }}
        tr:nth-child(even) {{ background-color: #f6f7f9; }}
        .title a {{ font-weight: bold; color: #0b5aa2; text-decoration: none; }}
        .title a:hover {{ text-decoration: underline; }}
        .topic {{ font-size: 16px; font-weight: bold; margin-top: 18px; }}
        .none {{ color: #777; font-style: italic; }}
        .note {{ font-size: 12px; color: #666; margin: 6px 0 10px; }}
        .query {{ margin: 6px 0 12px; font-size: 12px; color: #666; }}
        code {{ background: #f0f0f0; padding: 2px 4px; border-radius: 4px; }}
        .footer {{ margin-top: 18px; color: #666; font-size: 12px; }}
    </style>
    </head>
    <body>
      <h2>{DIGEST_TITLE}</h2>
      <div class="meta">
        Filter: <code>Published on {yesterday}</code> (yesterday only).<br/>
        Selection: <code>(anchor) AND (stats/method) AND (method-paper signals) AND NOT (clinical/discovery signals)</code>, then pick 1 paper per topic.
      </div>
    """

    for topic in ["Single-cell statistical methods", "Genomics statistical methods"]:
        payload = topic_payload.get(topic, {})
        article = payload.get("article")
        note = payload.get("note", "")
        query = payload.get("query", "")

        html += f'<div class="topic">{topic}</div>\n'
        if note:
            html += f'<div class="note">Selection note: <code>{note}</code></div>\n'
        if SHOW_QUERY and query:
            html += f'<div class="query">Query: <code>{query}</code></div>\n'

        html += """
        <table>
          <tr>
            <th style="width: 55%;">Title (click to read)</th>
            <th style="width: 25%;">First & Last Author</th>
            <th style="width: 20%;">Journal</th>
          </tr>
        """

        if article is None:
            html += f"""
              <tr>
                <td class="none" colspan="3">No paper found for this topic on {yesterday}.</td>
              </tr>
            """
        else:
            html += f"""
              <tr>
                <td class="title"><a href="{article['url']}">{article['title']}</a></td>
                <td>{article['authors']}</td>
                <td>{article['journal']}</td>
              </tr>
            """

        html += "</table>\n"

    html += """
      <div class="footer">
        Generated automatically by Yadi's PubMed Digest Agent (GitHub Actions + NCBI Entrez).
      </div>
    </body>
    </html>
    """
    return html


def send_email(html_content: str):
    if not EMAIL_SENDER or not EMAIL_PASSWORD or not EMAIL_RECEIVER:
        print("Error: EMAIL_SENDER / EMAIL_PASSWORD / EMAIL_RECEIVER not set (check GitHub Secrets).")
        return

    msg = MIMEMultipart()
    msg["From"] = EMAIL_SENDER
    msg["To"] = ", ".join(EMAIL_RECEIVER)
    msg["Subject"] = f"Yadi PubMed Digest (scRNA + Genomics methods) — {datetime.now().strftime('%Y-%m-%d')}"
    msg.attach(MIMEText(html_content, "html"))

    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.starttls()
        server.login(EMAIL_SENDER, EMAIL_PASSWORD)
        server.sendmail(EMAIL_SENDER, EMAIL_RECEIVER, msg.as_string())
        server.quit()
        print(f"Email successfully sent to {', '.join(EMAIL_RECEIVER)}")
    except Exception as e:
        print(f"Failed to send email: {e}")
        if "authentication" in str(e).lower():
            print("Authentication failed: ensure EMAIL_PASSWORD is a Google App Password.")


if __name__ == "__main__":
    if not EMAIL_PASSWORD or not ENTREZ_EMAIL:
        print("Agent could not run. Check EMAIL_PASSWORD and ENTREZ_EMAIL environment variables.")
    else:
        topic_payload = fetch_one_per_topic()
        html = format_html_email(topic_payload)
        send_email(html)
