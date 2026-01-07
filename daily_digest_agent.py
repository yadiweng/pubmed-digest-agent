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

# Convert the receiver string to a list of emails for smtplib.sendmail()
if EMAIL_RECEIVER_STR:
    EMAIL_RECEIVER = [e.strip() for e in EMAIL_RECEIVER_STR.split(",") if e.strip()]
else:
    EMAIL_RECEIVER = []

# =========================================================
# Topics: ONLY scRNA + genomics statistical methods
# Each topic will select at most 1 paper from yesterday
# =========================================================

# Topic-specific terms (domain anchor)
TOPICS = {
    "Single-cell statistical methods": [
        "single-cell RNA-seq", "scRNA-seq", "single-cell transcriptomics",
        "cell type annotation", "trajectory inference", "pseudotime",
        "batch correction", "data integration", "dimensionality reduction",
        "latent variable", "clustering"
    ],
    "Genomics statistical methods": [
        "genome-wide association", "GWAS", "statistical genetics",
        "fine-mapping", "polygenic risk score", "PRS",
        "eQTL", "genomic prediction", "rare variant",
        "high-dimensional", "genetic association"
    ]
}

# Method-oriented terms (to avoid pure biology/atlas papers)
# We keep it broad (Bayesian + high-dim + inference language).
METHOD_TERMS = [
    "statistical", "method", "methodology", "model", "modeling",
    "inference", "estimation", "likelihood",
    "bayesian", "posterior", "prior",
    "regression", "penalized", "lasso", "sparse",
    "variational", "expectation maximization", "EM",
    "probabilistic", "latent", "optimization"
]

# Max PMIDs per topic to fetch (we will pick 1)
CANDIDATE_RETMAX = 30


def _yesterday_pub_date_str() -> str:
    return (datetime.now() - timedelta(days=1)).strftime("%Y/%m/%d")


def _or_clause_titleab_or_mesh(terms: list[str]) -> str:
    """
    (("t1"[Title/Abstract] OR "t1"[MeSH Terms]) OR ( ... ))
    """
    clauses = [f'("{t}"[Title/Abstract] OR "{t}"[MeSH Terms])' for t in terms]
    return "(" + " OR ".join(clauses) + ")"


def _or_clause_titleab(terms: list[str]) -> str:
    """
    ("m1"[Title/Abstract] OR "m2"[Title/Abstract] OR ...)
    We use Title/Abstract for method terms to keep it interpretable and avoid over-broad MeSH.
    """
    clauses = [f'"{t}"[Title/Abstract]' for t in terms]
    return "(" + " OR ".join(clauses) + ")"


def _build_query(topic_terms: list[str], method_terms: list[str], pub_date: str) -> str:
    topic_clause = _or_clause_titleab_or_mesh(topic_terms)
    method_clause = _or_clause_titleab(method_terms)
    # Yesterday-only publication date filter
    return f"{topic_clause} AND {method_clause} AND {pub_date}[Date - Publication]"


def _search_pmids(query: str) -> list[str]:
    """
    Search PubMed and return PMIDs. Sort by pub date (within yesterday, still deterministic).
    """
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
                "pmid": pmid,
                "url": url,
            })
        except Exception as e:
            print(f"[WARN] Skipping article due to parsing error: {e}")
            continue

    return out


def fetch_one_per_topic() -> dict:
    """
    For each topic:
      - query = (topic terms) AND (method terms) AND (yesterday)
      - search PMIDs, fetch details
      - pick 1 representative paper (first returned, deterministic)
    """
    if not ENTREZ_EMAIL:
        print("Error: ENTREZ_EMAIL not set.")
        return {}

    Entrez.email = ENTREZ_EMAIL
    yesterday = _yesterday_pub_date_str()

    results = {}
    for topic, topic_terms in TOPICS.items():
        query = _build_query(topic_terms, METHOD_TERMS, yesterday)
        print(f"\n=== Topic: {topic} ===")
        print(f"Query: {query}")

        try:
            pmids = _search_pmids(query)
            if not pmids:
                print("No PMIDs found.")
                results[topic] = {"article": None, "query": query}
                continue

            articles = _fetch_article_details(pmids)
            if not articles:
                print("PMIDs found, but no parsable article details.")
                results[topic] = {"article": None, "query": query}
                continue

            results[topic] = {"article": articles[0], "query": query}
            print(f"Selected PMID {articles[0]['pmid']}")
        except Exception as e:
            print(f"[ERROR] Topic '{topic}' failed: {e}")
            results[topic] = {"article": None, "query": query}

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
        .query {{ margin: 6px 0 12px; font-size: 12px; color: #666; }}
        code {{ background: #f0f0f0; padding: 2px 4px; border-radius: 4px; }}
        .footer {{ margin-top: 18px; color: #666; font-size: 12px; }}
    </style>
    </head>
    <body>
      <h2>{DIGEST_TITLE}</h2>
      <div class="meta">
        Filter: <code>Published on {yesterday}</code> (yesterday only).<br/>
        Selection: For each topic, we apply <code>(topic terms) AND (method terms)</code>, then pick 1 representative paper.
      </div>
    """

    for topic in ["Single-cell statistical methods", "Genomics statistical methods"]:
        payload = topic_payload.get(topic, {})
        article = payload.get("article")
        query = payload.get("query", "")

        html += f'<div class="topic">{topic}</div>\n'
        if query:
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
