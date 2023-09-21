from pathlib import Path
from typing import Tuple, List, Union, Optional, Sequence

import abc
import bs4
import yattag
from bs4 import BeautifulSoup
from yattag import SimpleDoc


class HtmlBase(object):
    """
    Base class for HTML objects, contains common methods for classes that implement HTML objects.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        """
        Initializes an HTML base object.
        """
        self._doc = yattag.SimpleDoc()

    def get_tag(self, tag: str, attributes: List[Tuple[str, str]] = None) -> SimpleDoc.Tag:
        """
        Returns a tag with the given attributes.
        :param tag: Tag
        :param attributes: HTML attributes
        :return: Tag
        """
        if attributes:
            return self._doc.tag(tag, *attributes)
        else:
            return self._doc.tag(tag)

    def add_header(self, text: str, level: int, attributes: List[Tuple[str, str]] = None) -> None:
        """
        Adds a header.
        :param text: Header text
        :param level: Header level (1-6)
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        with self.get_tag('h{}'.format(str(level)), attributes):
            self.add_text(text)

    def add_paragraph(self, text: str, attributes: List[Tuple[str, str]] = None):
        """
        Adds a text paragraph
        :param text: Paragraph text
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        with self.get_tag('p', attributes):
            self.add_text(text)

    def add_line_break(self):
        """
        Adds a line break.
        :return: None
        """
        self._doc.stag('br')

    def to_html(self):
        """
        Converts the HTML object to HTML code.
        :return: HTML code
        """
        return self._doc.getvalue()

    def add_table(self, data: List[Sequence[Union[str, int, 'HtmlBase']]], column_names: List[str] = None,
                  table_attributes: List[Tuple[str, str]] = None) -> None:
        """
        Adds a table.
        :param data: Table data
        :param column_names: Column names
        :param table_attributes: Attributes of the table
        :return: None
        """
        with self.get_tag('table', table_attributes):
            if column_names is not None:
                self._add_table_header(column_names)
            for row in data:
                self._add_table_row(row)
        self.add_line_break()

    def _add_table_header(self, column_names: List[str]) -> None:
        """
        Adds a header to a table.
        :param column_names: Names of the columns
        :return: None
        """
        with self.get_tag('tr'):
            for name in column_names:
                with self.get_tag('th'):
                    self.add_text(name)

    def _add_table_row(self, data: Sequence[Union[str, int]], attributes: Optional[List[Tuple[str, str]]] = None) -> \
            None:
        """
        Adds a row to a table.
        :param data: Row data
        :param attributes: Attributes (optional)
        :return: None
        """
        with self.get_tag('tr', attributes):
            for value in data:
                if isinstance(value, HtmlBase):
                    self._doc.asis(value.to_html())
                else:
                    with self.get_tag('td'):
                        self.add_text(value)

    def add_link_to_file(self, link_text: str, relative_path: Path) -> None:
        """
        Adds a link to a file.
        :param link_text: Text of the link
        :param relative_path: Relative path to the file
        :return: None
        """
        if relative_path is None:
            self._doc.text(link_text)
            with self._doc.tag('em', klass='not-available'):
                self._doc.text('Not available')
        else:
            with self._doc.tag('a', href=str(relative_path)):
                self._doc.text(link_text)
        self.add_line_break()

    def add_error_message(self, message: str) -> None:
        """
        Adds an error message.
        :param message: Message text
        :return: None
        """
        self.add_alert(message, 'error')

    def add_warning_message(self, message: str) -> None:
        """
        Adds a warning message.
        :param message: Message text
        :return: None
        """
        self.add_alert(message, 'warning')

    def add_alert(self, message: str, type_: str) -> None:
        """
        Adds an alert to the report.
        Supported types: 'error', 'warning' and 'info'
        :param message: Message to display
        :param type_: Message type
        :return: None
        """
        with self._doc.tag('div', klass=f'alert-box {type_}'):
            with self._doc.tag('span'):
                self._doc.text(f'{type_}: ')
            self.add_text(message)

    def add_labeled_list(self, rows: List[Sequence[str]], ordered=False,
                         attributes: List[Tuple[str, str]] = None) -> None:
        """
        Adds a labeled list.
        :param rows: List of label-text pairs
        :param ordered: If True an ordered list is added
        :param attributes: Additional HTML arguments (key-value)
        :return: None
        """
        list_tag = 'ol' if ordered else 'ul'
        with self.get_tag(list_tag, attributes):
            for label, text in rows:
                with self._doc.tag('li'):
                    with self._doc.tag('b'):
                        self._doc.text('{}: '.format(label))
                    self.add_text(text)

    def add_html_object(self, input_object: 'HtmlBase') -> None:
        """
        Adds an HTML object.
        :param input_object: Input object
        :return: None
        """
        if not isinstance(input_object, HtmlBase):
            raise ValueError("{} is not an HTML object".format(input_object))
        self._doc.asis(input_object.to_html())

    def add_raw(self, html_code: str) -> None:
        """
        Adds raw HTML code, usage of this method is generally discouraged. Only when there is no other possibility.
        :param html_code: HTML code
        :return: None
        """
        self._doc.asis(html_code)

    def add_text(self, text: str) -> None:
        """
        Adds text to this HTML tag.
        :return: None
        """
        if text is None:
            raise ValueError("Cannot add None as text to the HTML object")
        html = BeautifulSoup(str(text), 'html.parser')
        for part in html.contents:
            if isinstance(part, bs4.element.Tag):
                self._doc.asis(f'<{part.name}>')
                self._doc.text(part.text)
                self._doc.asis(f'</{part.name}>')
            else:
                self._doc.text(str(part))

    def add_module_header(self, text: str, id_: str = None) -> None:
        """
        Adds a sub header.
        :param text: Header text
        :param id_: Html id
        :return: None
        """
        attributes = [('class', 'sub_header')]
        if id_ is not None:
            attributes.append(('id', id_))
        with self.get_tag('div', attributes):
            with self.get_tag('h2'):
                self.add_text(text)

    def add_horizontal_line(self) -> None:
        """
        Adds a horizontal line.
        :return: None
        """
        self._doc.stag('hr')
