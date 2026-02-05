import requests
from bs4 import BeautifulSoup
import csv
import os
from urllib.parse import urljoin
from datetime import datetime

URL = 'https://nipah-map.com/historical'
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'Data', 'Outbreak')
OUT_FILE = os.path.join(OUT_DIR, 'outbreaks_nipah.csv')

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36'
}


def text_or_empty(el):
    return el.get_text(strip=True) if el else ''


def parse_card(card, base_url=URL):
    # left column containing country, location, started, id, source
    left = None
    # try to find the left column by searching for anchor with /country/
    anchor = card.select_one('a[href^="/country/"]')
    country = ''
    country_href = ''
    if anchor:
        country = text_or_empty(anchor)
        country_href = urljoin(base_url, anchor.get('href'))

    # status may be in a span nearby
    status_el = card.select_one('div .inline-flex') or card.select_one('span')
    status = text_or_empty(status_el)

    # gather p tags under the left block (likely the first column)
    # find ancestor that contains both the anchor and p tags
    left = anchor
    while left is not None and left.name != 'div':
        left = left.parent

    # now scan nearby p tags under the same higher container
    location = ''
    started = ''
    nid = ''
    source = ''
    # prefer searching the card for 'Started:' / 'ID:' / 'Source:' patterns
    for p in card.find_all('p'):
        t = text_or_empty(p)
        if t.startswith('Started:'):
            started = t.replace('Started:', '').strip()
        elif t.startswith('ID:'):
            nid = t.replace('ID:', '').strip()
        elif t.startswith('Source:'):
            source = t.replace('Source:', '').strip()
        else:
            # first non-labeled p likely the location
            if not location:
                location = t

    # right column: cases, deaths, cfr
    cases = deaths = cfr = ''
    # find the stats container - look for children with numeric text
    # direct approach: find all divs with a child div that contains a numeric/text value
    stat_blocks = []
    # the stats are usually in a div with multiple child divs; find groups with 'Cases' label
    for block in card.find_all('div'):
        inner_div = block.find('div')
        if inner_div and inner_div.get_text(strip=True):
            # candidate: must also have a sibling label (like 'Cases')
            labels = [s.get_text(strip=True) for s in block.find_all('div')]
            if len(labels) >= 2:
                stat_blocks.append(labels)

    # stat_blocks may contain multiple repeats; pick the first reasonable 3-value group
    stats_vals = None
    for b in stat_blocks:
        # b example: ['2', 'Cases'] or flattened groups
        # try to detect groups of two where the second matches known labels
        if len(b) >= 2 and any(lbl in ['Cases', 'Deaths', 'CFR', 'CFR'] or 'Cases' in lbl for lbl in b[1:]):
            # accumulate numeric-like items from b
            nums = [x for x in b if x and not any(k in x for k in ['Cases', 'Deaths', 'CFR', '%'])]
            if nums:
                # fallback parsing below
                pass

    # Simpler targeted extraction: look for divs with class containing 'gap-4' (the stats container)
    stats_container = None
    for d in card.find_all('div'):
        cls = d.get('class')
        if cls and any('gap-4' in c for c in cls):
            stats_container = d
            break

    if stats_container:
        # direct children are three blocks (cases, deaths, cfr)
        children = [ch for ch in stats_container.find_all(recursive=False) if ch.name == 'div']
        vals = []
        for ch in children:
            v = ch.find('div')
            if v:
                vals.append(text_or_empty(v))
        if len(vals) >= 1:
            cases = vals[0]
        if len(vals) >= 2:
            deaths = vals[1]
        if len(vals) >= 3:
            cfr = vals[2]

    return {
        'country': country,
        'country_href': country_href,
        'status': status,
        'location': location,
        'started': started,
        'id': nid,
        'source': source,
        'cases': cases,
        'deaths': deaths,
        'cfr': cfr,
        'scraped_at': datetime.utcnow().isoformat()
    }


def scrape(url=URL):
    print('Requesting', url)
    r = requests.get(url, headers=HEADERS, timeout=15)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, 'html.parser')

    # find candidate card containers: look for divs that contain a country anchor and stats
    results = []
    anchors = soup.select('a[href^="/country/"]')
    seen = set()
    for a in anchors:
        # climb up to the outer card container: go up until a parent that contains both the anchor and a 'Cases' label
        node = a
        card = None
        for _ in range(6):
            parent = node.parent
            if parent is None:
                break
            # detect if parent contains a numeric stat div
            if parent.find('div', text=lambda t: t and 'Cases' in t) or parent.select_one('div[class*="gap-4"]'):
                card = parent
                break
            node = parent
        if card is None:
            # fallback: use the top-level ancestor 2 levels up
            card = a.find_parent('div')

        # use the card root (try to find the container with two child cols)
        root = card
        # avoid duplicates
        key = (a.get_text(strip=True), text_or_empty(root))
        if key in seen:
            continue
        seen.add(key)

        parsed = parse_card(root, base_url=url)
        results.append(parsed)

    return results


def save_csv(rows, out_file=OUT_FILE):
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    fieldnames = ['country', 'country_href', 'status', 'location', 'started', 'id', 'source', 'cases', 'deaths', 'cfr', 'scraped_at']
    with open(out_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    print('Saved', out_file)


def main():
    try:
        rows = scrape()
        if not rows:
            print('No cards found; the site may be JS-rendered or structure changed.')
        else:
            save_csv(rows)
    except Exception as e:
        print('Error:', e)


if __name__ == '__main__':
    main()
