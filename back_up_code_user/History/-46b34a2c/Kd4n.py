import requests
from bs4 import BeautifulSoup
from collections import Counter
import re

def get_article_links(url):
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    
    article_links = []
    for link in soup.find_all('a', href=True):
        if re.search(r'/doi/abs/\d+\.\d+\/PhysRevLett\.\d+\.\d+', link['href']):
            article_links.append('https://journals.aps.org' + link['href'])
    return article_links

def get_keywords_from_article(url):
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    
    keywords = []
    for keyword in soup.find_all('meta', attrs={'name': 'citation_keywords'}):
        keywords.extend([word.strip() for word in keyword['content'].split(',')])
    return keywords

def main():
    base_url = 'https://journals.aps.org/prl/issues/'
    recent_years = ['2022', '2023']

    keyword_counter = Counter()

    for year in recent_years:
        for month in range(1, 13):
            issue_url = f'{base_url}{year}/{month}'
            print(f'Fetching articles from {issue_url}')

            article_links = get_article_links(issue_url)
            for article_link in article_links:
                print(f'Getting keywords from {article_link}')
                keywords = get_keywords_from_article(article_link)
                keyword_counter.update(keywords)

    print('\nKeywords and their frequencies:')
    for keyword, frequency in keyword_counter.most_common():
        print(f'{keyword}: {frequency}')

if __name__ == '__main__':
    main()
