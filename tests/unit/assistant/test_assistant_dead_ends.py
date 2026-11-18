"""Regressions for answers that ended without a useful response.

Each test is one request that a user got no useful reply to although the
application had what it takes, with the fence unchanged: no processing setting is
touched, and a written answer is still checked against what was read.
"""
import json

import pytest

from app.schemas.assistant import AssistantRequest, Evidence, UiAction
from app.services import assistant_prose
from app.services import scientific_assistant as assistant
from app.services.assistant_methods import METHODS, WORKFLOW_GUIDE

from test_assistant_levels import _coverage_view
from test_assistant_tools_and_safety import _tool_turns

NO_RESULT = [Evidence(tool="current_plot", status="unavailable",
                      text="volcano_data|status=error|detail=Please select two different conditions.")]


@pytest.fixture
def advanced(assistant_rig):
    assistant_rig.profile["value"] = ("extended", 6, 2000)
    return assistant_rig


def _text_turns(monkeypatch, contents):
    """Model replies written as text instead of tool calls."""
    seen = []

    def turn(request, messages, tools, **kwargs):
        seen.append(json.loads(json.dumps(messages)))
        return "haiku", {"content": contents.pop(0)}
    monkeypatch.setattr(assistant, "native_tools_supported", lambda request: True)
    monkeypatch.setattr(assistant, "assistant_tool_turn", turn)
    return seen


# --- a page without a result -------------------------------------------------

def test_an_explanation_written_on_a_page_without_a_result_is_kept(advanced, monkeypatch):
    advanced.evidence = NO_RESULT
    _tool_turns(monkeypatch, [
        ("read_results", {}),
        ("answer", {"task": "current_plot", "text": "A moderated test shrinks each row's variance towards a shared estimate. "
                                                     "Nothing is calculated yet: select two different conditions first."})])
    result = advanced.run([], message="what does the moderated t-test do differently?", view=_coverage_view())
    assert result.answer.startswith("A moderated test shrinks")
    assert "No calculated result is available" not in result.answer


def test_that_explanation_is_still_checked_against_what_was_read(advanced, monkeypatch):
    advanced.evidence = NO_RESULT
    _tool_turns(monkeypatch, [
        ("read_results", {}),
        ("answer", {"task": "current_plot", "text": "Select two different conditions first. It will then show 4321 hits."})])
    result = advanced.run([], message="what does this page do?", view=_coverage_view())
    assert "Select two different conditions first." in result.answer and "4321" not in result.answer


def test_a_model_without_its_own_text_explains_from_the_method_and_the_guide(advanced):
    advanced.evidence = NO_RESULT
    plan = json.dumps({"task": "current_plot"})
    choice = json.dumps({"source": "result:qc.coverage", "query": ""})
    written = json.dumps({"answer": "Nothing is calculated yet. Use Download Plot beside the plot once it has rendered."})
    result = advanced.run([plan, choice, written], message="can this figure go into my thesis somehow?", view=_coverage_view())
    assert result.answer.startswith("Nothing is calculated yet.")
    sent = json.dumps(advanced.calls[-1])
    assert "No calculated result is available on the open page: Please select two different conditions" in sent
    assert "Download Plot" in sent  # the workflow guide travelled along


def test_a_level_that_picks_sentences_keeps_the_short_statement(assistant_rig):
    assistant_rig.evidence = NO_RESULT
    plan = json.dumps({"task": "current_plot"})
    choice = json.dumps({"source": "result:qc.coverage", "query": ""})
    result = assistant_rig.run([plan, choice], message="can this figure go into my thesis somehow?", view=_coverage_view())
    assert result.answer.startswith("No calculated result is available for this selection: Please select two different conditions")


def test_nothing_loaded_is_said_as_such(monkeypatch):
    monkeypatch.setattr(assistant, "get_all_current_datasets", lambda: {"protein": None, "phospho": None})
    request = AssistantRequest(message="what am I looking at?", view={"module": "qc.coverage"})
    found = assistant.scientific_evidence(request, assistant.AssistantPlan(task="current_plot"))
    assert found[0].status == "unavailable" and "No dataset is loaded yet" in found[0].text and "Data / Upload" in found[0].text
    assert assistant._unavailable_reason(found).startswith("No dataset is loaded yet")


def test_a_shortened_method_text_is_not_a_partly_read_result():
    note = assistant.partial_qualification(Evidence(tool="method:workflow.overview", status="partial", text="x"))
    assert "part of the result" in note  # the wording this rule keeps away from method texts
    source = open(assistant.__file__, encoding="utf-8").read()
    assert source.count('item.status == "partial" and not item.tool.startswith("method:")') == 2


# --- what the pages are --------------------------------------------------------

@pytest.mark.parametrize("page", ["data.upload", "data.annotation", "data.imputation", "data.restore", "settings.unknown"])
def test_pages_the_assistant_only_opens_are_described(page):
    assert page in METHODS and "assistant" in METHODS[page] or page == "data.restore"
    assert assistant._page_summary(page)


def test_the_notebook_export_is_described_as_it_is_built():
    for text in (WORKFLOW_GUIDE, METHODS["summary.report"]):
        assert "captured input" not in text and "no data matrix" in text
    assert "Ship data matrices" in WORKFLOW_GUIDE and "Ship full imputed matrix" in WORKFLOW_GUIDE
    assert "Save State" in WORKFLOW_GUIDE and "Ctrl/Cmd+S" in METHODS["data.restore"]
    assert METHODS["data.distribution"].startswith("Distribution summary and QQ plot")


@pytest.mark.parametrize("message,page", [
    ("where can I compare two samples against each other?", "comparison.pearson"),
    ("where do I upload a peptide report?", "data.upload"),
    ("where can I switch on that the notebook ships the data matrices?", "settings.unknown"),
])
def test_where_questions_find_those_pages(assistant_rig, message, page):
    pages = {"stats.anova": "Statistical Analysis / One-way ANOVA", "summary.report": "Summary / Report",
             "comparison.pearson": "Comparison / Pearson Correlation", "data.upload": "Data / Upload", "settings.unknown": "Settings"}
    result = assistant_rig.run([], message=message, view={"module": "qc.coverage", "destinations": pages})
    assert result.answer.split("\n")[1].startswith(f"- {pages[page]}") and assistant_rig.calls == []


# --- kinase activity ---------------------------------------------------------

def test_a_tool_planning_model_cannot_promise_kinase_activity(advanced, monkeypatch):
    _tool_turns(monkeypatch, [("answer", {"task": "current_plot", "module": "phospho.ksea",
                                         "text": "Once the comparison is set I can list the most active kinases."})])
    advanced.evidence = None
    monkeypatch.setattr(assistant, "scientific_evidence", ORIGINAL_EVIDENCE)
    result = advanced.run([], message="which kinases are active in GRK?", view={"module": "stats.volcano", "workflow": "phosphosite"})
    assert "inputs, not kinase-activity results" in result.answer and "upload both result files" in result.answer
    assert "most active kinases" not in result.answer


ORIGINAL_EVIDENCE = assistant.scientific_evidence


@pytest.mark.parametrize("message,asks", [
    ("which kinases are active in GRK?", True), ("can I see kinase activity somewhere?", True),
    ("what is upstream of these sites, which kinase?", True), ("open the kinase library", False),
    ("is the AKT1 site S473 changed?", False)])
def test_what_counts_as_a_kinase_activity_question(message, asks):
    assert bool(assistant._KINASE_ACTIVITY.search(message)) is asks


# --- replies in the wrong form -------------------------------------------------

def test_answer_arguments_written_as_text_are_that_answer(advanced, monkeypatch):
    _text_turns(monkeypatch, [json.dumps({"task": "answer", "module": "data.restore",
                                          "text": "Use Save State on Data / Restore."})])
    # Written without a read, it is answered from the page's own description (one writing call).
    written = json.dumps({"answer": "Use Save State on Data / Restore."})
    result = advanced.run([written], message="Wie speichere ich meine Sitzung?", view=_coverage_view())
    assert advanced.plans[0].module == "data.restore" and advanced.plans[0].task in {"current_plot", "help"}
    assert "Save State" in result.answer and "Save State (also Ctrl/Cmd+S)" in json.dumps(advanced.calls[-1])


def test_an_empty_change_list_written_as_text_gets_the_one_reminder(advanced, monkeypatch):
    seen = _text_turns(monkeypatch, [json.dumps({"task": "ui", "text": "", "modules": []}),
                                     json.dumps({"task": "ui", "actions": [{"type": "navigate", "target": "qc.coverage"}]})])
    advanced.run([], message="bring me to the coverage plot", view={**_coverage_view(), "module": "qc.cv"})
    assert len(seen) == 2 and seen[1][-1]["content"] == assistant._USE_A_TOOL


# --- processing decisions stay with the user, and the user is told where --------

def test_a_refusal_in_the_models_own_words_is_not_lost(advanced, monkeypatch):
    _tool_turns(monkeypatch, [("answer", {"task": "current_plot", "text": "I can't remove a sample; that changes the data."})])
    result = advanced.run([], message="remove sample GRK_5, it looks like an outlier", view=_coverage_view())
    assert result.needsClarification and result.actions == []
    assert "Data / Annotation" in result.answer and "No supported action" not in result.answer


def test_an_unknown_answer_style_does_not_lose_the_answer(assistant_rig):
    plan = json.dumps({"task": "current_plot", "responseStyle": "direct"})
    choice = json.dumps({"source": "result:qc.coverage", "query": ""})
    assistant_rig.evidence = NO_RESULT
    result = assistant_rig.run([plan, choice], message="can this figure go into my thesis somehow?", view=_coverage_view())
    assert result.answer.startswith("No calculated result is available for this selection")


def test_a_long_refusal_keeps_its_whole_sentences(advanced, monkeypatch):
    text = "I haven't removed GRK_5. " + "Dropping a sample changes the data behind the plot, so I can't do it from here. " * 6
    _tool_turns(monkeypatch, [("answer", {"task": "current_plot", "text": text})])
    result = advanced.run([], message="remove sample GRK_5, it looks like an outlier", view=_coverage_view())
    assert result.needsClarification and result.actions == []
    assert result.answer.startswith("I haven't removed GRK_5. Dropping a sample") and "Data / Annotation" in result.answer
    assert "so I can't do it from here. Which samples belong" in result.answer


@pytest.mark.parametrize("message,answer,added", [
    ("filter out proteins with less than 5 valid values", "I can't apply a filter from this plot.", "Data / Annotation"),
    ("set the localization cutoff to 0.75", "I can't change it; use the Peptide Collapse tool.", "field shown with the phosphosite plots"),
    ("change the imputation to the numerical method", "I can't change imputation. Set it on Data / Imputation.", ""),
    ("is the median higher in GRK?", "I can't tell from this plot.", ""),
    ("normalize with VSN", "Done.", ""),
])
def test_a_refused_processing_change_names_the_users_page(message, answer, added):
    result = assistant._pipeline_place_added(message, answer)
    assert (result != answer) is bool(added) and added in result and result.startswith(answer)


def test_a_model_question_naming_another_page_gives_way_to_the_applications():
    wrong = "Should I open Data / Imputation so you can set the filter?"
    result = assistant._with_pipeline_place("filter out proteins with less than 5 valid values", wrong)
    assert "Data / Imputation" not in result and "Data / Annotation" in result and "never changes" in result
    kept = assistant._with_pipeline_place("impute the missing values", "Imputation is yours to run on Data / Imputation.")
    assert kept == "Imputation is yours to run on Data / Imputation."


# --- names, levels and looks ---------------------------------------------------

def _entity_view(module="singleProtein.boxplot", workflow="protein", **mode):
    return {"module": module, "workflow": workflow, "mode": {"kind": "phospho" if workflow == "phosphosite" else workflow, **mode},
            "destinations": {module: "Single Protein"}, "controls": [
                {"id": "e", "label": "Protein / Site", "type": "select", "value": "A_HUMAN", "options": ["A_HUMAN", "B_HUMAN"],
                 "optionCount": 7450, "permission": "confirm"}]}


def test_a_gene_name_the_list_does_not_show_is_resolved_from_the_table(advanced, monkeypatch):
    monkeypatch.setattr(assistant, "feature_candidates", lambda kind, parent="": ["G3P_HUMAN"] if parent == "GAPDH" else [])
    plan = json.dumps({"task": "ui", "actions": [{"type": "search_options", "target": "e", "value": "GAPDH"}]})
    result = advanced.run([plan], message="now the boxplot for a single protein, GAPDH", view=_entity_view())
    assert [(a.type, a.value) for a in result.actions] == [("set_control", "G3P_HUMAN")]


def test_a_synonym_the_model_searched_gives_way_to_the_name_the_user_wrote(advanced, monkeypatch):
    monkeypatch.setattr(assistant, "feature_candidates", lambda kind, parent="": ["ARBK1_HUMAN"] if parent == "GRK2" else [])
    plan = json.dumps({"task": "ui", "actions": [{"type": "search_options", "target": "e", "value": "ADRBK1"}]})
    result = advanced.run([plan], message="show GRK2 as a boxplot", view=_entity_view())
    assert [(a.type, a.value) for a in result.actions] == [("set_control", "ARBK1_HUMAN")]


def test_a_name_the_user_did_not_give_is_not_resolved(advanced, monkeypatch):
    monkeypatch.setattr(assistant, "feature_candidates", lambda kind, parent="": ["G3P_HUMAN"])
    plan = json.dumps({"task": "ui", "actions": [{"type": "search_options", "target": "e", "value": "GAPDH"}]})
    result = advanced.run([plan], message="show me something interesting", view=_entity_view())
    assert [(a.type, a.value) for a in result.actions] == [("search_options", "GAPDH")]


def test_its_phosphosites_selects_the_parent_protein_of_the_named_one(advanced, monkeypatch):
    monkeypatch.setattr(assistant, "feature_candidates",
                        lambda kind, parent="": ["P31749~AKT1_S124_M1", "P31749~AKT1_S129_M1"] if parent == "AKT1" else [])
    plan = json.dumps({"task": "ui", "actions": [{"type": "search_options", "target": "e", "value": "AKT1"}]})
    result = advanced.run([plan], message="show its phosphosites as a heatmap", view=_entity_view("singleProtein.heatmap", "phosphosite"),
                          messages=[{"role": "user", "content": "show me AKT1"}, {"role": "assistant", "content": "Showing AKT1_HUMAN."},
                                    {"role": "user", "content": "show its phosphosites as a heatmap"}])
    assert [(a.type, a.value) for a in result.actions] == [("set_control", "P31749")]


def test_the_dataset_level_goes_first_and_alone(advanced):
    view = {"module": "stats.volcano", "workflow": "protein", "destinations": {"stats.volcano": "Volcano"}, "controls": [
        {"id": "lvl", "label": "Dataset level", "type": "select", "value": "protein", "options": ["protein", "phospho"], "permission": "view"},
        {"id": "c1", "label": "Condition 1", "type": "select", "value": "WT", "options": ["WT", "GRK"], "permission": "confirm"},
        {"id": "c2", "label": "Condition 2", "type": "select", "value": "WT", "options": ["WT", "GRK"], "permission": "confirm"}]}
    plan = json.dumps({"task": "ui", "actions": [{"type": "set_control", "target": "c1", "value": "GRK"},
                                                 {"type": "set_control", "target": "lvl", "value": "phospho"}]})
    result = advanced.run([plan], message="volcano of the phosphosites, GRK vs WT", view=view)
    assert [(a.target, a.value) for a in result.actions] == [("lvl", "phospho")] and not result.finished


@pytest.mark.parametrize("message", ["the sample labels are too small", "the x tick labels are tiny", "Die Beschriftung ist zu klein"])
def test_too_small_asks_for_a_size(message):
    view = _coverage_view()
    request = AssistantRequest(message=message, view=view)
    control = next(control for control in request.view.controls if control.id == "xt")
    assert assistant._scope_problem(request, UiAction(type="set_control", target="xt", value="12"), control) is None


def test_save_a_look_and_use_it_everywhere_are_two_steps(assistant_rig):
    message = "save this look as paper and use it for all plots"
    saved = assistant_rig.run([], message=message, view=_coverage_view())
    assert [(a.type, a.value) for a in saved.actions] == [("save_style", "paper")] and not saved.finished
    applied = assistant_rig.run([], message=message, view=_coverage_view(), savedStyles=["paper"],
                                receipts=['Saved this plot\'s look as the style "paper".'])
    assert json.loads(applied.actions[0].value) == {"name": "paper", "scope": "all"} and applied.finished
    assert assistant_rig.calls == []


def test_localization_counts_and_their_shares_are_confirmed():
    evidence = ["project_localization_counts|sites_total=11580|sites_at_or_above_0.50=7153|sites_at_or_above_0.75=5174"]
    text = ("Of the 11,580 sites, 7,153 reach a probability of 0.50 (about 62%). "
            "At 0.75 there are 5,174 sites (about 45%). About 80% are above 0.5.")
    checked, removed = assistant_prose.verify(text, evidence, "how well are the sites localized")
    assert "7,153" in checked and "about 45%" in checked and "80%" not in checked and removed == 1


def test_a_settings_search_on_a_page_without_controls_returns_what_the_page_holds(advanced, monkeypatch):
    seen = _tool_turns(monkeypatch, [("find_settings", {"word": "Rscript"}),
                                    ("answer", {"task": "help", "module": "settings.unknown",
                                                "text": "R Compatibility on this page sets the R installation used for VSN and limma."})])
    written = json.dumps({"answer": "R Compatibility on this page sets the R installation used for VSN and limma."})
    result = advanced.run([written], message="where do I set which R is used?",
                          view={"module": "settings.unknown", "destinations": {"settings.unknown": "Settings"}, "controls": []})
    assert "R Compatibility" in seen[1]["messages"][-1]["content"]
    assert result.answer.startswith("R Compatibility on this page")


# --- the read-only viewer (2026-11-13) -----------------------------------------

@pytest.fixture
def viewer(monkeypatch):
    monkeypatch.setenv("COPYLOT_VIEWER_MODE", "1")


def test_the_viewer_describes_its_data_pages_as_read_only(viewer):
    from app.services.assistant_methods import VIEWER_NOTE, method_evidence

    assert method_evidence("data.upload").text == VIEWER_NOTE and "cannot be changed here" in VIEWER_NOTE
    assert method_evidence("workflow.overview").text.startswith(VIEWER_NOTE)
    assert method_evidence("qc.cv").text == METHODS["qc.cv"]  # a plot is explained as in the analysis app


def test_outside_the_viewer_nothing_of_it_is_said():
    from app.services.assistant_methods import VIEWER_NOTE, method_evidence

    assert method_evidence("data.upload").text == METHODS["data.upload"]
    assert VIEWER_NOTE not in method_evidence("workflow.overview").text


def test_every_read_in_the_viewer_says_what_a_viewer_cannot_do(viewer, monkeypatch):
    monkeypatch.setattr(assistant, "_scientific_evidence",
                        lambda *args: [Evidence(tool="current_plot", status="available", text="volcano|up=7|down=5")])
    found = assistant.scientific_evidence(AssistantRequest(message="can I download the notebook?", view={"module": "stats.volcano"}),
                                          assistant.AssistantPlan(task="current_plot"))
    assert found[-1].tool == "viewer" and "Notebook exports and save states are not offered" in found[-1].text


def test_a_refused_processing_change_in_the_viewer_names_no_page(viewer):
    refused = assistant._with_pipeline_place("normalize the data please", "Should I open Data / Annotation for you?")
    assert "Data / Annotation" not in refused and "read-only viewer" in refused
    added = assistant._pipeline_place_added("normalize the data please",
                                            "I can't change normalization. Check Data / Annotation for it.")
    assert added.startswith("I can't change normalization.") and "Data / Annotation" not in added and "read-only viewer" in added
    assert assistant._pipeline_place_added("normalize the data please", refused) == refused  # said once


def test_a_named_comparison_is_still_to_come_after_the_level(advanced):
    view = {"module": "stats.volcano", "workflow": "protein", "destinations": {"stats.volcano": "Volcano"}, "controls": [
        {"id": "lvl", "label": "Dataset level", "type": "select", "value": "protein", "options": ["protein", "phospho"], "permission": "view"},
        {"id": "c1", "label": "Condition 1", "type": "select", "value": "C0", "options": ["C0", "C5", "dQ0"], "permission": "confirm"},
        {"id": "c2", "label": "Condition 2", "type": "select", "value": "C0", "options": ["C0", "C5", "dQ0"], "permission": "confirm"}]}
    plan = json.dumps({"task": "ui", "actions": [{"type": "set_control", "target": "lvl", "value": "phospho"}]})
    named = advanced.run([plan], message="volcano of the phosphosites, C5 vs C0", view=view)
    assert [a.target for a in named.actions] == ["lvl"] and not named.finished
    alone = advanced.run([plan], message="switch this volcano to the phospho data", view=view)
    assert [a.target for a in alone.actions] == ["lvl"] and alone.finished
