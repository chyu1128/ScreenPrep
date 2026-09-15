from __future__ import annotations

import csv
import os
import tkinter as tk
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, FilterCatalog
    from rdkit.Chem.SaltRemover import SaltRemover
    from PIL import ImageTk
except ImportError as exc:  # pragma: no cover - handled at launch
    raise SystemExit(
        "RDKit is required. Install it with: conda install -c conda-forge rdkit pandas"
    ) from exc


REQUIRED_COLUMNS = ("ID", "SMILES", "VALUE")


@dataclass
class Record:
    compound_id: str
    original_smiles: str
    value: str
    cleaned_smiles: str = ""
    outcome: str = "Pending"
    detail: str = ""
    initially_excluded: bool = False


class ScreenPrep(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("ScreenPrep · Drug Screening Preprocessor")
        self.geometry("1180x760")
        self.minsize(940, 620)
        self.configure(bg="#f5f8fa")

        self.records: list[Record] = []
        self.source_path: Path | None = None
        self.output_records: list[Record] = []
        self.missing_records: list[Record] = []
        self.duplicate_records: list[Record] = []
        self.pains_records: list[Record] = []
        self.retained_record_by_duplicate: dict[int, Record] = {}
        self.table_item_by_record: dict[int, str] = {}
        self.record_by_table_item: dict[str, Record] = {}
        self.missing_index = 0
        self.processed = False
        self.pains_var = tk.BooleanVar(value=True)
        self.path_var = tk.StringVar(value="No file selected")
        self.summary_var = tk.StringVar(value="")
        self._setup_style()
        self._build_ui()

    def _setup_style(self) -> None:
        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure("App.TFrame", background="#f5f8fa")
        style.configure("Card.TFrame", background="white")
        style.configure("Title.TLabel", background="#f5f8fa", foreground="#102a43", font=("Segoe UI", 23, "bold"))
        style.configure("Subtitle.TLabel", background="#f5f8fa", foreground="#627d98", font=("Segoe UI", 10))
        style.configure("CardTitle.TLabel", background="white", foreground="#102a43", font=("Segoe UI", 11, "bold"))
        style.configure("CardText.TLabel", background="white", foreground="#627d98", font=("Segoe UI", 9))
        style.configure("Primary.TButton", background="#0b7285", foreground="white", font=("Segoe UI", 10, "bold"), padding=(16, 10))
        style.map("Primary.TButton", background=[("active", "#095c6b")])
        style.configure("Secondary.TButton", background="#e7f5f7", foreground="#0b7285", font=("Segoe UI", 9, "bold"), padding=(12, 8))
        style.configure("Keep.TButton", background="#d3f9d8", foreground="#1b5e20", font=("Segoe UI", 9, "bold"), padding=(12, 8))
        style.map("Keep.TButton", background=[("active", "#b2f2bb")])
        style.configure("Exclude.TButton", background="#ffe3e3", foreground="#b42318", font=("Segoe UI", 9, "bold"), padding=(12, 8))
        style.map("Exclude.TButton", background=[("active", "#ffc9c9")])
        style.configure("Treeview", rowheight=31, font=("Segoe UI", 9), background="white", fieldbackground="white", foreground="#243b53")
        style.configure("Treeview.Heading", font=("Segoe UI", 9, "bold"), background="#edf2f7", foreground="#486581", relief="flat")
        style.map("Treeview", background=[("selected", "#d9f0f4")], foreground=[("selected", "#102a43")])

    def _build_ui(self) -> None:
        shell = ttk.Frame(self, style="App.TFrame", padding=(32, 25))
        shell.pack(fill="both", expand=True)

        header = ttk.Frame(shell, style="App.TFrame")
        header.pack(fill="x")
        ttk.Label(header, text="ScreenPrep", style="Title.TLabel").pack(anchor="w")
        ttk.Label(header, text="Prepare clean, screen-ready structures from a raw assay CSV", style="Subtitle.TLabel").pack(anchor="w", pady=(2, 20))

        controls = ttk.Frame(shell, style="Card.TFrame", padding=20)
        controls.pack(fill="x")
        controls.columnconfigure(1, weight=1)
        ttk.Label(controls, text="1  INPUT DATA", style="CardTitle.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(controls, textvariable=self.path_var, style="CardText.TLabel").grid(row=0, column=1, sticky="w", padx=(18, 10))
        ttk.Button(controls, text="Choose CSV", style="Secondary.TButton", command=self.load_csv).grid(row=0, column=2, sticky="e")
        ttk.Label(controls, text="Expected columns:  ID, SMILES, VALUE", style="CardText.TLabel").grid(row=1, column=0, columnspan=3, sticky="w", pady=(10, 0))

        line = ttk.Separator(controls, orient="horizontal")
        line.grid(row=2, column=0, columnspan=3, sticky="ew", pady=16)
        ttk.Label(controls, text="2  PREPROCESSING", style="CardTitle.TLabel").grid(row=3, column=0, sticky="w")
        ttk.Checkbutton(controls, text="Remove PAINS hits", variable=self.pains_var, command=self._refresh_filter_note).grid(row=3, column=1, sticky="w", padx=(18, 10))
        self.filter_note = ttk.Label(controls, text="Enabled · RDKit PAINS catalog", style="CardText.TLabel")
        self.filter_note.grid(row=3, column=2, sticky="e")
        ttk.Label(controls, text="Duplicate structures are always removed after salt removal.", style="CardText.TLabel").grid(row=4, column=0, columnspan=3, sticky="w", pady=(10, 0))

        actions = ttk.Frame(shell, style="App.TFrame")
        actions.pack(fill="x", pady=(18, 14))
        ttk.Button(actions, text="Run preprocessing", style="Primary.TButton", command=self.process).pack(side="left")
        self.export_button = ttk.Button(actions, text="Export cleaned CSV", style="Secondary.TButton", command=self.export_csv, state="disabled")
        self.export_button.pack(side="left", padx=10)
        self.issue_actions = ttk.Frame(actions, style="App.TFrame")
        self.issue_actions.pack(side="left")
        summary_box = ttk.Frame(actions, style="App.TFrame")
        summary_box.pack(side="right", fill="x", expand=True, padx=(16, 0))
        ttk.Label(summary_box, textvariable=self.summary_var, style="Subtitle.TLabel", anchor="e").pack(fill="x", expand=True)

        table_card = ttk.Frame(shell, style="Card.TFrame", padding=(18, 15))
        table_card.pack(fill="both", expand=True)
        ttk.Label(table_card, text="Processing preview", style="CardTitle.TLabel").pack(anchor="w")
        ttk.Label(table_card, text="Missing values, invalid structures, PAINS hits, and duplicate structures are excluded from export.", style="CardText.TLabel").pack(anchor="w", pady=(3, 12))
        self.review_panel = ttk.Frame(table_card, style="Card.TFrame")
        self.review_actions = ttk.Frame(self.review_panel, style="Card.TFrame")
        for index, width in enumerate((95, 230, 230, 95, 130, 240)):
            self.review_actions.columnconfigure(index, minsize=width)
        self.structure_button = ttk.Button(self.review_actions, text="Structure\n ", width=12, style="Secondary.TButton", command=self.show_structure, state="disabled")
        self.show_retained_button = ttk.Button(self.review_actions, text="Show retain\nrecord", width=12, style="Secondary.TButton", command=self.show_retained_record, state="disabled")
        self.keep_button = ttk.Button(self.review_actions, text="Keep\n ", width=12, style="Keep.TButton", command=self.keep_selected, state="disabled")
        self.exclude_button = ttk.Button(self.review_actions, text="Exclude\n ", width=12, style="Exclude.TButton", command=self.exclude_selected, state="disabled")
        self.structure_button.grid(row=0, column=2, sticky="w", padx=4)
        self.show_retained_button.grid(row=0, column=3, sticky="w", padx=4)
        self.keep_button.grid(row=0, column=4, sticky="w", padx=4)
        self.exclude_button.grid(row=0, column=5, sticky="w", padx=4)
        self.review_log_box = ttk.Frame(self.review_panel, style="Card.TFrame")
        ttk.Label(self.review_log_box, text="Manual review log", style="CardText.TLabel").pack(anchor="w", pady=(0, 4))
        log_holder = ttk.Frame(self.review_log_box, style="Card.TFrame")
        log_holder.pack(fill="x")
        self.review_log = tk.Text(
            log_holder,
            height=3,
            wrap="none",
            state="disabled",
            relief="solid",
            borderwidth=1,
            bg="#f8fafc",
            fg="#486581",
            font=("Consolas", 9),
        )
        log_scroll = ttk.Scrollbar(log_holder, orient="vertical", command=self.review_log.yview)
        self.review_log.configure(yscrollcommand=log_scroll.set)
        self.review_log.pack(side="left", fill="x", expand=True)
        log_scroll.pack(side="right", fill="y")
        self.review_panel.pack(fill="x")
        columns = ("id", "input", "cleaned", "value", "outcome", "detail")
        self.table = ttk.Treeview(table_card, columns=columns, show="headings", selectmode="browse")
        labels = {"id": "ID", "input": "Input SMILES", "cleaned": "Cleaned SMILES", "value": "VALUE", "outcome": "Outcome", "detail": "Notes"}
        widths = {"id": 95, "input": 230, "cleaned": 230, "value": 95, "outcome": 130, "detail": 240}
        for key in columns:
            self.table.heading(key, text=labels[key])
            self.table.column(key, width=widths[key], minwidth=70, anchor="w")
        scroll_y = ttk.Scrollbar(table_card, orient="vertical", command=self.table.yview)
        scroll_x = ttk.Scrollbar(table_card, orient="horizontal", command=self.table.xview)
        self.table.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
        self.table.bind("<<TreeviewSelect>>", self._on_table_selection)
        self.table.pack(side="left", fill="both", expand=True)
        scroll_y.pack(side="right", fill="y")
        scroll_x.pack(side="bottom", fill="x")

        footer = ttk.Label(shell, text="RDKit-based preprocessing · Salt removal uses RDKit's standard salt definitions", style="Subtitle.TLabel")
        footer.pack(anchor="w", pady=(13, 0))

    def _refresh_filter_note(self) -> None:
        self.filter_note.configure(text="Enabled · RDKit PAINS catalog" if self.pains_var.get() else "Disabled")

    def load_csv(self) -> None:
        filename = filedialog.askopenfilename(title="Choose screening CSV", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not filename:
            return
        try:
            with open(filename, "r", newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                if not reader.fieldnames or not set(REQUIRED_COLUMNS).issubset(reader.fieldnames):
                    raise ValueError("CSV must contain the columns: ID, SMILES, VALUE")
                self.records = [
                    Record(
                        (row.get("ID", "") or "").strip(),
                        (row.get("SMILES", "") or "").strip(),
                        (row.get("VALUE", "") or "").strip(),
                    )
                    for row in reader
                ]
            self.source_path = Path(filename)
            self.path_var.set(f"{self.source_path.name}  ·  {len(self.records):,} records")
            self.summary_var.set("")
            self.output_records = []
            self.missing_records = []
            self.duplicate_records = []
            self.pains_records = []
            self.retained_record_by_duplicate = {}
            self.processed = False
            self._clear_review_log()
            self.export_button.configure(state="disabled")
            self._update_issue_actions()
            self._update_review_actions()
            self._render_rows(self.records)
        except (OSError, ValueError, csv.Error) as exc:
            messagebox.showerror("Cannot read CSV", str(exc))

    def process(self) -> None:
        if not self.records:
            messagebox.showinfo("Choose an input", "Select a CSV file before preprocessing.")
            return
        running_dialog = self._show_running_dialog()
        remover = SaltRemover()
        pains = FilterCatalog.FilterCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        kept: list[Record] = []
        seen_structures: dict[str, Record] = {}
        self.missing_records = []
        self.duplicate_records = []
        self.pains_records = []
        self.retained_record_by_duplicate = {}
        self.missing_index = 0
        self._clear_review_log()
        missing = invalid = salts = pains_removed = duplicates_removed = 0
        for record in self.records:
            record.cleaned_smiles, record.outcome, record.detail = "", "", ""
            record.initially_excluded = False
            empty_fields = [name for name, value in zip(REQUIRED_COLUMNS, (record.compound_id, record.original_smiles, record.value)) if not value]
            if empty_fields:
                record.outcome, record.detail = "Excluded", f"Missing value: {', '.join(empty_fields)}"
                missing += 1
                self.missing_records.append(record)
                continue
            mol = Chem.MolFromSmiles(record.original_smiles)
            if mol is None:
                record.outcome, record.detail = "Excluded", "Invalid SMILES"
                invalid += 1
                continue
            desalted = remover.StripMol(mol, dontRemoveEverything=True)
            record.cleaned_smiles = Chem.MolToSmiles(desalted, canonical=True)
            original_canonical = Chem.MolToSmiles(mol, canonical=True)
            was_desalted = original_canonical != record.cleaned_smiles
            pains_match = pains.GetFirstMatch(desalted) if self.pains_var.get() else None
            if pains_match:
                record.outcome, record.detail = "Excluded", f"PAINS alert: {pains_match.GetDescription()}"
                pains_removed += 1
                self.pains_records.append(record)
                continue
            if record.cleaned_smiles in seen_structures:
                record.outcome = "Excluded"
                retained_record = seen_structures[record.cleaned_smiles]
                record.detail = f"Duplicate of {retained_record.compound_id}"
                duplicates_removed += 1
                self.duplicate_records.append(record)
                self.retained_record_by_duplicate[id(record)] = retained_record
                continue
            seen_structures[record.cleaned_smiles] = record
            record.outcome = "Ready"
            record.detail = "Salt removed" if was_desalted else "Structure validated"
            salts += int(was_desalted)
            kept.append(record)
        for record in self.records:
            record.initially_excluded = record.outcome == "Excluded"
        self.output_records = kept
        self._render_rows(self.records)
        self._update_issue_actions()
        self.export_button.configure(state="normal" if kept else "disabled")
        self.processed = True
        self._update_review_actions()
        self.summary_var.set(f"{len(kept):,} ready  ·  {missing} missing  ·  {invalid} invalid  ·  {duplicates_removed} duplicates  ·  {salts} desalted  ·  {pains_removed} PAINS")
        if running_dialog.winfo_exists():
            running_dialog.destroy()
        messagebox.showinfo("Preprocessing complete", "Preprocessing is complete. Review the results and export the cleaned CSV when ready.", parent=self)

    def _show_running_dialog(self) -> tk.Toplevel:
        """Display a lightweight, non-dismissible indicator while RDKit processes records."""
        dialog = tk.Toplevel(self)
        dialog.title("Preprocessing")
        dialog.transient(self)
        dialog.resizable(False, False)
        dialog.protocol("WM_DELETE_WINDOW", lambda: None)
        ttk.Frame(dialog, padding=(28, 22)).pack(fill="both", expand=True)
        content = dialog.winfo_children()[0]
        ttk.Label(content, text="Preprocessing is running…", font=("Segoe UI", 11, "bold")).pack()
        ttk.Label(content, text="Validating structures, removing salts, and applying selected filters.", style="Subtitle.TLabel", justify="center").pack(pady=(8, 0))
        dialog.update_idletasks()
        dialog.update()
        return dialog

    def _render_rows(self, rows: list[Record]) -> None:
        self.table.delete(*self.table.get_children())
        self.table_item_by_record = {}
        self.record_by_table_item = {}
        for record in rows:
            tag = "excluded" if record.outcome == "Excluded" else "ready"
            item = self.table.insert("", "end", values=(record.compound_id, record.original_smiles, record.cleaned_smiles, record.value, record.outcome, record.detail), tags=(tag,))
            self.table_item_by_record[id(record)] = item
            self.record_by_table_item[item] = record
        self.table.tag_configure("excluded", foreground="#b42318")
        self.table.tag_configure("ready", foreground="#087f5b")

    def _update_review_actions(self) -> None:
        if self.processed:
            self.review_actions.pack(fill="x", pady=(0, 12))
            self.review_log_box.pack(fill="x", pady=(0, 12))
        else:
            self.review_actions.pack_forget()
            self.review_log_box.pack_forget()
        for button in (self.structure_button, self.show_retained_button, self.keep_button, self.exclude_button):
            button.configure(state="disabled")

    def _on_table_selection(self, _event: object = None) -> None:
        if not self.processed:
            return
        state = "normal" if self.table.selection() else "disabled"
        for button in (self.structure_button, self.show_retained_button, self.keep_button, self.exclude_button):
            button.configure(state=state)

    def _clear_review_log(self) -> None:
        self.review_log.configure(state="normal")
        self.review_log.delete("1.0", "end")
        self.review_log.configure(state="disabled")

    def _append_review_log(self, action: str, record: Record) -> None:
        self.review_log.configure(state="normal")
        self.review_log.insert("end", f"{action} {record.compound_id or '(blank ID)'}\n")
        self.review_log.see("end")
        self.review_log.configure(state="disabled")

    def _selected_record(self) -> Record | None:
        selection = self.table.selection()
        return self.record_by_table_item.get(selection[0]) if selection else None

    def show_retained_record(self) -> None:
        record = self._selected_record()
        if not record:
            messagebox.showerror("Select a duplicate", "Select a duplicate row before showing its retained record.")
            return
        retained = self.retained_record_by_duplicate.get(id(record))
        if not retained:
            messagebox.showerror("No retained record", "The selected row is not an excluded duplicate.")
            return
        item = self.table_item_by_record.get(id(retained))
        if item:
            self.table.focus(item)
            self.table.selection_set(item)
            self.table.see(item)
            self.table.xview_moveto(0)
            self.table.focus_set()

    def show_structure(self) -> None:
        """Render the selected record's cleaned structure with RDKit."""
        record = self._selected_record()
        if not record:
            messagebox.showerror("Select a record", "Select a record before viewing its molecular structure.")
            return
        smiles = record.cleaned_smiles or record.original_smiles
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            messagebox.showerror("No structure available", "The selected row does not contain a valid SMILES structure.")
            return
        dialog = tk.Toplevel(self)
        dialog.title(f"Molecular structure · {record.compound_id or '(blank ID)'}")
        dialog.configure(bg="white")
        image = ImageTk.PhotoImage(Draw.MolToImage(mol, size=(560, 380)))
        ttk.Label(dialog, text=record.compound_id or "(blank ID)", font=("Segoe UI", 11, "bold")).pack(padx=20, pady=(18, 6))
        molecule = ttk.Label(dialog, image=image, background="white")
        molecule.image = image
        molecule.pack(padx=20, pady=(0, 12))
        ttk.Label(dialog, text=smiles, wraplength=540, justify="center", background="white", foreground="#486581").pack(padx=20, pady=(0, 18))

    def keep_selected(self) -> None:
        record = self._selected_record()
        if not record:
            return
        if record.outcome != "Excluded":
            messagebox.showerror("Record already kept", "The selected record is already included in the export.")
            return
        if not record.cleaned_smiles:
            messagebox.showerror("Cannot keep record", "Rows with missing values or invalid SMILES cannot be retained for export.")
            return
        if not messagebox.askyesno("Keep compound", f"Keep {record.compound_id or '(blank ID)'} in the export?", parent=self):
            return
        record.outcome, record.detail = "Ready", "Kept manually"
        self.missing_records = [item for item in self.missing_records if item is not record]
        self.duplicate_records = [item for item in self.duplicate_records if item is not record]
        self.pains_records = [item for item in self.pains_records if item is not record]
        self._append_review_log("Kept", record)
        self._refresh_after_manual_change(record)

    def exclude_selected(self) -> None:
        record = self._selected_record()
        if not record:
            return
        if record.initially_excluded:
            messagebox.showerror("Originally excluded", "This record was originally excluded by preprocessing and cannot be manually excluded again.")
            return
        if record.outcome == "Excluded":
            messagebox.showerror("Record already excluded", "The selected record is already excluded from export.")
            return
        if not messagebox.askyesno("Exclude compound", f"Exclude {record.compound_id or '(blank ID)'} from the export?", parent=self):
            return
        record.outcome, record.detail = "Excluded", "Excluded manually"
        self._append_review_log("Excluded", record)
        self._refresh_after_manual_change(record)

    def _refresh_after_manual_change(self, selected: Record) -> None:
        self.output_records = [record for record in self.records if record.outcome == "Ready"]
        self._render_rows(self.records)
        self._update_issue_actions()
        self.export_button.configure(state="normal" if self.output_records else "disabled")
        item = self.table_item_by_record.get(id(selected))
        if item:
            self.table.focus(item)
            self.table.selection_set(item)
            self.table.see(item)
        self._on_table_selection()

    def _update_issue_actions(self) -> None:
        """Show issue controls only when the current run found those issues."""
        for child in self.issue_actions.winfo_children():
            child.destroy()
        if self.missing_records:
            ttk.Button(
                self.issue_actions,
                text=f"Show missing ({len(self.missing_records)})",
                style="Secondary.TButton",
                command=self.show_missing,
            ).pack(side="left", padx=(0, 8))
        if self.duplicate_records:
            ttk.Button(
                self.issue_actions,
                text=f"Show duplicates ({len(self.duplicate_records)})",
                style="Secondary.TButton",
                command=self.show_duplicates,
            ).pack(side="left")
        if self.pains_records:
            ttk.Button(
                self.issue_actions,
                text=f"Show PAINS ({len(self.pains_records)})",
                style="Secondary.TButton",
                command=self.show_pains,
            ).pack(side="left", padx=(8, 0))

    def show_missing(self) -> None:
        """Select a missing record and bring its empty field's column into view."""
        if not self.missing_records:
            return
        record = self.missing_records[self.missing_index % len(self.missing_records)]
        self.missing_index += 1
        item = self.table_item_by_record.get(id(record))
        if not item:
            return
        self.table.focus(item)
        self.table.selection_set(item)
        self.table.see(item)
        values = (record.compound_id, record.original_smiles, record.value)
        missing_column = next(index for index, value in enumerate(values) if not value)
        # The displayed input fields are ID, Input SMILES, Cleaned SMILES, then VALUE.
        self.table.xview_moveto((0.0, 0.10, 0.45)[missing_column])
        self.table.focus_set()

    def show_duplicates(self) -> None:
        """Present each excluded duplicate with the original retained compound."""
        if not self.duplicate_records:
            return
        dialog = tk.Toplevel(self)
        dialog.title("Duplicate structures")
        dialog.geometry("820x360")
        dialog.minsize(650, 260)
        dialog.configure(bg="#f5f8fa")
        panel = ttk.Frame(dialog, style="App.TFrame", padding=22)
        panel.pack(fill="both", expand=True)
        ttk.Label(panel, text="Duplicate structures", style="Title.TLabel").pack(anchor="w")
        ttk.Label(panel, text="Compared using canonical SMILES after salt removal. The first record is retained.", style="Subtitle.TLabel").pack(anchor="w", pady=(2, 3))
        table_holder = ttk.Frame(panel, style="App.TFrame")
        table_holder.pack(fill="both", expand=True)
        columns = ("duplicate", "retained", "structure")
        table = ttk.Treeview(table_holder, columns=columns, show="headings", selectmode="browse")
        for name, text, width in (
            ("duplicate", "Excluded duplicate", 175),
            ("retained", "Retained record", 160),
            ("structure", "Canonical cleaned SMILES", 430),
        ):
            table.heading(name, text=text)
            table.column(name, width=width, minwidth=100, anchor="w")
        for record in self.duplicate_records:
            retained_record = self.retained_record_by_duplicate[id(record)]
            table.insert("", "end", values=(record.compound_id or "(blank ID)", retained_record.compound_id or "(blank ID)", record.cleaned_smiles))
        scrollbar = ttk.Scrollbar(table_holder, orient="vertical", command=table.yview)
        table.configure(yscrollcommand=scrollbar.set)
        table.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        ttk.Button(panel, text="Close", style="Secondary.TButton", command=dialog.destroy).pack(anchor="e", pady=(14, 0))

    def show_pains(self) -> None:
        """Present PAINS-filtered structures from the current preprocessing run."""
        if not self.pains_records:
            return
        dialog = tk.Toplevel(self)
        dialog.title("PAINS structural alerts")
        dialog.geometry("820x360")
        dialog.minsize(650, 260)
        dialog.configure(bg="#f5f8fa")
        panel = ttk.Frame(dialog, style="App.TFrame", padding=22)
        panel.pack(fill="both", expand=True)
        ttk.Label(panel, text="PAINS structural alerts", style="Title.TLabel").pack(anchor="w")
        ttk.Label(panel, text="These records were excluded by the optional RDKit PAINS catalog filter.", style="Subtitle.TLabel").pack(anchor="w", pady=(2, 16))
        holder = ttk.Frame(panel, style="App.TFrame")
        holder.pack(fill="both", expand=True)
        columns = ("id", "structure", "alert")
        table = ttk.Treeview(holder, columns=columns, show="headings", selectmode="browse")
        for name, text, width in (("id", "ID", 150), ("structure", "Cleaned SMILES", 390), ("alert", "PAINS alert", 240)):
            table.heading(name, text=text)
            table.column(name, width=width, minwidth=100, anchor="w")
        for record in self.pains_records:
            table.insert("", "end", values=(record.compound_id or "(blank ID)", record.cleaned_smiles, record.detail.removeprefix("PAINS alert: ")))
        scrollbar = ttk.Scrollbar(holder, orient="vertical", command=table.yview)
        table.configure(yscrollcommand=scrollbar.set)
        table.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        ttk.Button(panel, text="Close", style="Secondary.TButton", command=dialog.destroy).pack(anchor="e", pady=(14, 0))

    def export_csv(self) -> None:
        if not self.output_records:
            return
        default_name = f"{self.source_path.stem if self.source_path else 'screening'}_cleaned.csv"
        filename = filedialog.asksaveasfilename(title="Export cleaned screening data", defaultextension=".csv", initialfile=default_name, filetypes=[("CSV files", "*.csv")])
        if not filename:
            return
        try:
            with open(filename, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS)
                writer.writeheader()
                writer.writerows({"ID": row.compound_id, "SMILES": row.cleaned_smiles, "VALUE": row.value} for row in self.output_records)
            messagebox.showinfo("Export complete", f"Saved {len(self.output_records):,} screening-ready records.\n\n{os.path.basename(filename)}")
        except OSError as exc:
            messagebox.showerror("Export failed", str(exc))


if __name__ == "__main__":
    ScreenPrep().mainloop()
