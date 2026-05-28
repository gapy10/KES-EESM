import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QComboBox, QPushButton, QLabel


def run_gui():
    app = QApplication(sys.argv)
    selected = {"index": None, "value": None}

    class DropdownApp(QWidget):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("Select an Option")
            self.setGeometry(100, 100, 300, 120)

            layout = QVBoxLayout()
            self.setLayout(layout)
            self.raw_options = ["Inducirana napetost", "Karakteristika prostega teka", "Statorsko napajanje", "Navor", "Ld", "Lq"]
            self.display_options = [f"{i + 1}. {opt}" for i, opt in enumerate(self.raw_options)]

            self.combo = QComboBox()
            self.combo.addItems(self.display_options)
            layout.addWidget(self.combo)

            self.button = QPushButton("Confirm")
            self.button.clicked.connect(self.confirm_selection)
            self.button.setFixedWidth(100)
            layout.addWidget(self.button)

        def confirm_selection(self):
            idx = self.combo.currentIndex()
            selected["index"] = idx + 1  # 1-based
            selected["value"] = self.raw_options[idx]
            self.close()  # Closes the window

    window = DropdownApp()
    window.show()
    app.exec_()

    return selected["index"], selected["value"]

def show_value_window(naloga, value):
    app = QApplication(sys.argv)

    # Create main window
    window = QWidget()
    window.setWindowTitle(naloga)
    window.resize(320, 160)

    # Create label with the variable's value
    label = QLabel(value)

    # Create close button
    close_button = QPushButton("Zapri")
    close_button.clicked.connect(window.close)
    close_button.setFixedWidth(100)

    # Layout setup
    layout = QVBoxLayout()
    layout.addWidget(label)
    layout.addWidget(close_button)
    window.setLayout(layout)

    # Show window
    window.show()
    sys.exit(app.exec_())