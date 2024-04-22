import json
from importlib.resources import files
from typing import Dict, Any

from pacu.app.report.htmlbase import HtmlBase


class HtmlCitation(HtmlBase):
    """
    This class is used to format citations in HTML.
    :param citation_data: Citation data (parsed from RIS)
    :return: None
    """

    def __init__(self, citation_data: Dict[str, Any]) -> None:
        """
        Initializes the citation.
        :param citation_data: Citation data
        """
        super().__init__()
        self._citation_data = citation_data
        self._authors = '; '.join(citation_data['first_authors'])
        self._pub_year = citation_data['publication_year'].split('/')[0]
        self._title = citation_data['primary_title']
        citation_processor = {
            'jour': self._process_journal_citation,
            'book': self._process_book_citation}
        try:
            import pprint
            pprint.pprint(citation_data['type_of_reference'].lower())
            citation_processor[citation_data['type_of_reference'].lower()]()
        except KeyError as err:
            raise KeyError(
                f"Citations of type {citation_data['type_of_reference']} are not supported! (KeyError '{err}')")

    def _process_journal_citation(self) -> None:
        """
        Creates the html tags when the citation is for a journal article.
        :return: None
        """
        with self.get_tag('div', attributes=[('class', 'citations')]):
            with self.get_tag('p'):
                # Volume + (number optional)
                journal_parts = [
                    self._citation_data['alternate_title3'], f", {self._citation_data.get('volume', 'n/a')}"]
                if 'number' in self._citation_data:
                    journal_parts.append(f" ({self._citation_data['number']})")
                journal = ''.join(journal_parts)

                # Citation text
                self.add_text(f"{self._authors} ({self._pub_year}). {self._title}. In <i>{journal}</i>. ")

                # Citation DOI / link
                with self.get_tag('a', [('href', f"https://dx.doi.org/{self._citation_data['doi']}")]):
                    self.add_text(f"DOI: {self._citation_data['doi']}")

    def _process_book_citation(self) -> None:
        """
        Creates the html tags when the citation is for a book.
        :return: None
        """
        with self.get_tag('div', attributes=[('class', 'citations')]):
            with self.get_tag('p'):
                edition = f", {self._citation_data['edition']}" if 'edition' in self._citation_data else ''
                publisher = self._citation_data['publisher']
                self.add_text(f"{self._authors} ({self._pub_year}). {self._title}{edition}. {publisher}")

    @staticmethod
    def parse_from_json(json_basename: str) -> 'HtmlCitation':
        """
        Parses a HTML citation from a JSON file.
        :param json_basename: Basename for the JSON file
        :return: Citation
        """
        json_citation = files('pacu').joinpath(f'resources/citations/{json_basename}.json')
        with json_citation.open(encoding='utf-8') as handle:
            data = json.load(handle)
        return HtmlCitation(data)
