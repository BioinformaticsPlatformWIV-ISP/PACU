import bs4
from bs4 import BeautifulSoup
from yattag import Doc

from pacu.app.report.htmlbase import HtmlBase


class HtmlElement(HtmlBase):
    """
    This class contains a custom HTML element.
    """

    def __init__(self, tag, text=None, attributes=None):
        """
        Initializes an HTML element.
        :param tag: Tag
        :param attributes: HTML attributes
        :param text: Text
        """
        super(HtmlElement, self).__init__()
        self._tag_name = tag
        self._attributes = attributes
        self._tag_text = text

    @property
    def text(self):
        """
        Returns the HTML element as plain text.
        :return: Text
        """
        return self._tag_text

    def _has_nested_content(self):
        """
        Returns true if this HTML element has nested content.
        :return: True if there is nested content
        """
        return len(self._doc.getvalue()) != 0

    # noinspection PyArgumentList
    def to_html(self):
        """
        Converts this element to HTML code. A novel Doc() instance is created in order to nest the content of this
        elements Doc() inside the tag associated with this HtmlElement.
        :return: HTML code
        """
        if self._attributes is None:
            self._attributes = []
        doc, tag, text = Doc().tagtext()
        if not self._has_nested_content() and self._tag_text is None:
            # If there is no content / text in the tag, a self closing tag (stag) is created
            doc.stag(self._tag_name, *self._attributes)
        else:
            # Otherwise the content / text is added inside of the tag
            with tag(self._tag_name, *self._attributes):
                if self._tag_text is not None:
                    html = BeautifulSoup(str(self._tag_text), 'html.parser')
                    for part in html.contents:
                        if isinstance(part, bs4.element.Tag):
                            doc.asis(f'<{part.name}>')
                            doc.text(part.text)
                            doc.asis(f'</{part.name}>')
                        else:
                            # noinspection PyTypeChecker
                            doc.text(part)
                if self._has_nested_content():
                    doc.asis(self._doc.getvalue())
        return doc.getvalue()
