import logging
from typing import Optional, List

import shutil

from pathlib import Path

from pacu.app.report.htmlbase import HtmlBase
from pacu.app.report.htmlelement import HtmlElement


class HtmlReport(HtmlBase):
    """
    This class represents an HTML report.
    """

    def __init__(self, filename: Path, output_dir: Path = None, include_js: Optional[List[Path]] = None):
        """
        Initializes the report.
        :param filename: Filename
        :param output_dir: Output directory
        :param include_js: (Optional) List of Javascript files that are included in the report
        """
        super(HtmlReport, self).__init__()
        self._filename = filename
        self._output_dir = output_dir
        self._include_js = include_js

    @property
    def output_dir(self) -> Path:
        """
        Returns the output directory.
        :return: Output directory
        """
        return self._output_dir

    def save(self) -> None:
        """
        Saves the report.
        :return: None
        """
        logging.info(f"Saving report '{self._filename.name}'")
        if self._filename is None:
            raise ValueError("Report with filename 'None' cannot be saved")
        with self._filename.open('w', encoding='utf-8') as handle:
            self.add_raw('<!DOCTYPE HTML>')
            handle.write(self.to_html())
        if self._include_js:
            if self._output_dir is None:
                raise ValueError("Cannot enable Javascript when there is no output directory")
            for path in self._include_js:
                shutil.copyfile(path, self._output_dir / path.name)

    def initialize(self, title: str, css_style: Path = None) -> None:
        """
        Initializes an HTML report.
        :param title: Report title
        :param css_style: CSS style
        :return: None
        """
        logging.info("Initializing report")
        with self._doc.tag('head'):
            with self._doc.tag('title'):
                self._doc.text(title)
            self._doc.stag('meta', charset='UTF-8')
            if self._include_js:
                for path in self._include_js:
                    with self.get_tag('script', [('type', 'text/javascript'), ('src', path.name)]):
                        self._doc.text('')
            if css_style is not None:
                with css_style.open() as css, self._doc.tag('style', type="text/css"):
                    self._doc.asis(css.read())

    def add_pipeline_header(self, pipeline_name: str) -> None:
        """
        Adds the pipeline header to the report.
        :param pipeline_name: Name of the pipeline
        :return: None
        """
        if self._output_dir is None:
            raise ValueError("Can't add the pipeline header without an output directory")
        with self.get_tag('div', [('class', 'header')]):
            with self.get_tag('div', [('id', 'header_title')]):
                self.add_text(f"{pipeline_name} report")

    def to_html(self) -> str:
        """
        Returns the report HTML code.
        :return: HTML code.
        """
        parent = HtmlElement('html')
        parent.add_raw(self._doc.getvalue())
        return parent.to_html()
