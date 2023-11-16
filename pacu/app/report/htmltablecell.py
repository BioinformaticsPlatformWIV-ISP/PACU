from pacu.app.report.htmlelement import HtmlElement


class HtmlTableCell(HtmlElement):
    """
    This class represents a cell in an HTML table.
    """

    def __init__(self, text, color=None, attributes=None, link=None):
        """
        Initializes a table cell.
        :param text: Text
        :param color: Color
        :param attributes: Attributes
        :param link: Link to add to the cell text
        """
        self._as_text = text
        if color is not None:
            color_attribute = [('class', color)]
            if attributes is not None:
                attributes += color_attribute
            else:
                attributes = color_attribute
        if link is not None:
            super().__init__('td', None, attributes)
            self.add_html_object(HtmlElement('a', text, [('href', link)]))
        else:
            super().__init__('td', text, attributes)

    @property
    def text(self):
        """
        Returns the text belonging to this tag.
        :return: Text
        """
        return self._as_text
