import logging
import shutil
from pathlib import Path
from typing import List, Union, Tuple, Optional, Sequence

from pacu.app.report.htmlbase import HtmlBase
from pacu.app.report.htmlelement import HtmlElement
from pacu.app.report.htmltablecell import HtmlTableCell


class HtmlReportSection(HtmlElement):
    """
    This class can be used to create a section in the HTML report.
    """

    def __init__(self, title: Optional[str], level: int = 2, subtitle: Optional[str] = None) -> None:
        """
        Initializes a report section.
        :param title: Section title
        :param level: Header level
        :param subtitle: Subtitle to put next to the header
        """
        super(HtmlReportSection, self).__init__('div', attributes=[('class', 'report_section')])
        self._title = title
        if (title is not None) and (subtitle is not None):
            self.add_header_with_subtitle(title, level, subtitle)
        elif title is not None:
            self.add_header(title, level, subtitle)
        self._files = []

    def add_header_with_subtitle(self, header: str, level: int, subtitle: str) -> None:
        """
        Adds a header with a subtitle.
        :param header: Header
        :param level: Header level
        :param subtitle: Subtitle
        :return: None
        """
        with self.get_tag('h{}'.format(str(level))):
            self.add_text(header + ' ')
            with self.get_tag('small', [('class', 'header-subtitle')]):
                self.add_text('({})'.format(subtitle))

    @property
    def files(self) -> List[Path]:
        """
        Returns the files that were added to this report.
        :return: Files
        """
        return self._files

    def copy_files(self, output_directory: Path) -> None:
        """
        Exports the files belonging to this section to the given directory.
        :param output_directory: Output directory
        :return: None
        """
        logging.info("Exporting report section files")
        if output_directory is None or not output_directory.is_dir():
            raise ValueError(f'Invalid output directory: {output_directory}')
        for file_path, relative_path in self.files:
            logging.info(f'Copying file: {file_path} -> {relative_path}')
            if not file_path.is_file():
                raise ValueError(f"Cannot add file (does not exist) '{file_path}'")
            relative_dir = output_directory / relative_path.parent
            if not relative_dir.is_dir():
                relative_dir.mkdir(parents=True)
            shutil.copy(file_path, output_directory / relative_path)

    def add_file(self, input_file: Path, relative_path: Path) -> Path:
        """
        Adds the file to the report.
        Returns the relative path parameter so it is easier to use in list comprehensions.
        :param input_file: Input file
        :param relative_path: path where the file will be saved relative to the report output directory
        :return: Relative path
        """
        logging.info("Adding file to report section: {}".format(relative_path))
        self._files.append((input_file, relative_path,))
        return relative_path

    def add_table(self, data: List[Sequence[Union[str, int, 'HtmlBase']]], column_names: List[str] = None,
                  table_attributes: List[Tuple[str, str]] = None, msg_if_empty: str = 'None found') -> None:
        """
        Adds a table to the report section.
        :param data: Table data
        :param column_names: Column names
        :param table_attributes: Table attributes
        :param msg_if_empty: If set, this message is added when the data is empty
        :return: None
        """
        if (msg_if_empty is not None) and (len(data) == 0) and (column_names is not None):
            data.append([HtmlTableCell(msg_if_empty, attributes=[('colspan', len(column_names))])])
        super().add_table(data, column_names, table_attributes)

    def __repr__(self) -> str:
        """
        Returns the internal representation.
        :return: String representation
        """
        return f"HtmlReportSection(title='{self._title}')"
