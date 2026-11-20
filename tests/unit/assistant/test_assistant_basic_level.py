"""What the application does for a Basic-level model only (qwen3.5:2b, 2026-11-13).

A small local model is guided more closely: fewer evidence sources to choose
from, and three requests the application settles itself. The Standard and
Advanced levels keep their own handling, which each test's counterpart checks.
"""
import json

import pytest

from app.schemas.assistant import AssistantRequest
from app.services import assistant_safety
from app.services import scientific_assistant as assistant

DESTINATIONS = {
    "qc.coverage": "QC Pipeline / Coverage Plot", "qc.pca": "QC Pipeline / PCA Plot",
    "data.distribution": "Data / Distribution", "data.restore": "Data / Save State",
    "data.translator": "Data / ID Translator", "stats.volcano": "Statistical Analysis / Volcano Plot",
    "phospho.ksea": "Phospho-specific / Kinase Activity",
}


def _view(module="qc.coverage", **extra):
    return {"module": module, "destinations": DESTINATIONS, **extra}


def _levels(rig, level):
    rig.profile["value"] = {"compact": ("compact", 2, 600), "standard": ("standard", 4, 1000),
                            "extended": ("extended", 6, 2000)}[level]


def _offered(rig):
    return next(record["data"]["sources"] for record in rig.records if "sources" in record["data"])


# --- the evidence choice -----------------------------------------------------

def test_a_basic_model_chooses_among_the_open_page_and_the_named_pages(assistant_rig):
    plan = json.dumps({"task": "help"})
    choice = json.dumps({"source": "method:data.restore", "query": ""})
    assistant_rig.run([plan, choice, "[]"], message="how do I save my session?", view=_view())
    offered = _offered(assistant_rig)
    assert "method:data.restore" in offered and "method:qc.coverage" in offered and "current_result" in offered
    assert "method:data.translator" not in offered and "method:stats.volcano" not in offered
    assert {"unsupported", "entity_lookup", "external_activity", "project_state", "method:workflow.overview"} <= set(offered)
    assert len(offered) < 20


@pytest.mark.parametrize("level", ["standard", "extended"])
def test_the_larger_levels_keep_every_source(assistant_rig, monkeypatch, level):
    _levels(assistant_rig, level)
    monkeypatch.setattr(assistant, "native_tools_supported", lambda request: False)
    plan = json.dumps({"task": "help"})
    choice = json.dumps({"source": "method:data.restore", "query": ""})
    assistant_rig.run([plan, choice, "[]", "[]"], message="how do I save my session?", view=_view())
    offered = _offered(assistant_rig)
    assert "method:data.translator" in offered and "method:stats.volcano" in offered and len(offered) > 60


# --- an order to process the data --------------------------------------------

@pytest.mark.parametrize("message, place", [
    ("normalize the data with the median please", "Data / Annotation"),
    ("can you impute the missing values", "Data / Imputation"),
    ("please remove the outlier sample WT_3", "Data / Annotation"),
])
def test_a_processing_order_is_declined_by_the_application(assistant_rig, message, place):
    result = assistant_rig.run([], message=message, view=_view())
    assert assistant_rig.calls == [] and not result.actions
    assert result.answer.startswith(assistant_safety.REFUSAL) and place in result.answer
    assert result.answer.endswith("Nothing was changed.")


@pytest.mark.parametrize("message", [
    "how was the data normalized?", "which normalization should I use?", "remove the legend",
    "normalize the colour scale", "filter the table for GRK2", "normalize rows by z-score",
])
def test_questions_and_display_requests_are_not_orders(message):
    assert not assistant._PIPELINE_ORDER.search(message)


def test_the_larger_levels_word_that_refusal_themselves(assistant_rig, monkeypatch):
    _levels(assistant_rig, "standard")
    request = AssistantRequest(message="normalize the data with the median please", model="synthetic", view=_view())
    assert assistant._basic_level(request) is False
    monkeypatch.setattr(assistant, "_basic_pipeline_response",
                        lambda *args: pytest.fail("only the Basic level is answered by the application"))
    plan = json.dumps({"task": "clarify", "clarification": "pipeline"})
    result = assistant_rig.run([plan], message=request.message, view=_view())
    assert assistant_rig.calls and "Data / Annotation" in result.answer


# --- kinase activity ---------------------------------------------------------

def test_a_basic_model_cannot_answer_kinase_activity_with_another_result(assistant_rig):
    plan = json.dumps({"task": "current_plot"})
    assistant_rig.run([plan, "{}", "{}", "{}"], message="which kinases are active in GRK?", view=_view("stats.volcano"))
    assert assistant_rig.plans and assistant_rig.plans[0].module == "external_activity"
    assert not any("sources" in call for call in assistant_rig.calls)


# --- a page that was asked for -----------------------------------------------

def test_an_open_request_that_came_back_as_nothing_to_do_opens_the_named_page(assistant_rig):
    view = _view(controls=[{"id": "known", "label": "Known", "type": "number", "value": "45"}])
    plan = json.dumps({"task": "ui", "actions": [{"type": "set_control", "target": "known", "value": "45.0"}]})
    result = assistant_rig.run([plan], message="open the QQ plot", view=view)
    assert [(action.type, action.target) for action in result.actions] == [("navigate", "data.distribution")]
    assert result.finished and not result.answer


def test_the_same_reply_at_the_standard_level_is_reported_as_before(assistant_rig):
    _levels(assistant_rig, "standard")
    view = _view(controls=[{"id": "known", "label": "Known", "type": "number", "value": "45"}])
    plan = json.dumps({"task": "ui", "actions": [{"type": "set_control", "target": "known", "value": "45.0"}]})
    result = assistant_rig.run([plan], message="open the QQ plot", view=view)
    assert not result.actions and result.answer == "These proposed view changes are already in place."


def test_a_question_on_the_open_page_is_not_sent_to_a_page_sharing_a_word(assistant_rig):
    request = AssistantRequest(message="which test should I use, my samples are paired", model="synthetic",
                               view=_view("stats.volcano"))
    offered = assistant._basic_sources(request, {f"method:{key}": "" for key in DESTINATIONS} | {"current_result": ""})
    assert set(offered) == {"method:stats.volcano", "current_result"}


def test_kinase_evidence_is_only_for_a_kinase_question(assistant_rig):
    plan = json.dumps({"task": "current_plot", "query": "sites"})
    choice = json.dumps({"source": "external_activity", "query": "how many phosphosites did we find?"})
    assistant_rig.run([plan, choice, "{}", "{}", "{}"], message="how many phosphosites did we find? and where is PCA",
                      view=_view("stats.volcano"))
    assert assistant_rig.plans[0].task == "current_plot" and assistant_rig.plans[0].module in {"", "stats.volcano"}


# --- a protein or comparison named without a verb ------------------------------

def _boxplot_view(value="A1"):
    return {"module": "singleProtein.boxplot", "workflow": "protein", "destinations": DESTINATIONS,
            "controls": [{"id": "p", "label": "Protein", "type": "select", "value": value, "permission": "confirm",
                          "options": ["A1", "B2"], "optionCount": 7450}]}


def test_the_named_protein_is_selected_on_the_page_that_was_opened(assistant_rig, monkeypatch):
    monkeypatch.setattr(assistant, "feature_candidates", lambda kind, parent="": ["G3P_HUMAN"] if parent == "GAPDH" else [])
    result = assistant_rig.run([], message="now the boxplot for a single protein, GAPDH", view=_boxplot_view(),
                               receipts=["Opened Single Protein / Protein Boxplot."])
    assert assistant_rig.calls == []
    assert [(action.type, action.target, action.value) for action in result.actions] == [("set_control", "p", "G3P_HUMAN")]


@pytest.mark.parametrize("found", [[], ["G3P_HUMAN", "A1"]])
def test_a_name_that_is_not_exactly_one_uploaded_protein_is_left_to_the_model(assistant_rig, monkeypatch, found):
    monkeypatch.setattr(assistant, "feature_candidates", lambda kind, parent="": found)
    request = AssistantRequest(message="boxplot, GAPDH", model="synthetic", view=_boxplot_view())
    assert assistant._basic_feature_response(request, 1) is None


def test_a_question_about_a_protein_selects_nothing(assistant_rig, monkeypatch):
    monkeypatch.setattr(assistant, "feature_candidates", lambda kind, parent="": ["G3P_HUMAN"])
    request = AssistantRequest(message="is GAPDH significant?", model="synthetic", view=_boxplot_view())
    assert assistant._basic_feature_response(request, 1) is None


def test_a_page_named_without_open_is_opened(assistant_rig):
    result = assistant_rig.run([], message="volcano plot please", view=_view())
    assert assistant_rig.calls == []
    assert [(action.type, action.target) for action in result.actions] == [("navigate", "stats.volcano")]


@pytest.mark.parametrize("message", ["PCA looks weird", "what does the volcano plot show?", "why is the PCA like this"])
def test_a_remark_or_question_about_a_page_does_not_open_it(message):
    request = AssistantRequest(message=message, model="synthetic", view=_view())
    assert assistant._basic_destination_response(request, 1) is None


# --- text that is too small ---------------------------------------------------

def _sized_view():
    def size(id_, label):
        return {"id": id_, "label": label, "type": "number", "value": "10", "appearance": True, "permission": "view"}
    return _view(controls=[size("t", "Title / Size"), size("xt", "X axis / X tick label size"), size("lg", "Legend / Text size")])


def test_too_small_without_a_named_size_lists_the_sizes(assistant_rig):
    result = assistant_rig.run([], message="the sample labels are too small", view=_sized_view())
    assert assistant_rig.calls == [] and not result.actions and result.needsClarification
    assert "title, x tick label, legend text" in result.answer and result.answer.endswith("Nothing was changed.")


def test_too_small_with_the_size_named_makes_it_bigger(assistant_rig):
    result = assistant_rig.run([], message="the title is too small", view=_sized_view())
    assert assistant_rig.calls == []
    assert [(action.target, action.value) for action in result.actions] == [("t", "12.0")]


def test_too_small_does_not_pick_a_setting_by_one_shared_word(assistant_rig):
    view = _sized_view()
    view["controls"].append({"id": "vl", "label": "Value labels / Label size", "type": "number", "value": "8",
                             "appearance": True, "permission": "view"})
    result = assistant_rig.run([], message="the sample labels are too small", view=view)
    assert not result.actions and result.needsClarification and "value labels label" in result.answer


# --- downloads and replicate quality ------------------------------------------

@pytest.mark.parametrize("message", ["how do I get this plot out for my paper?", "can I export this figure as SVG?",
                                     "how do I download the table?"])
def test_how_to_download_is_the_guides_sentence(assistant_rig, message):
    result = assistant_rig.run([], message=message, view=_view("stats.volcano"))
    assert assistant_rig.calls == [] and not result.actions
    assert "Download Plot" in result.answer and "Download table (.xlsx)" in result.answer and "Export notebook" in result.answer
    assert "which saves an Excel workbook" in result.answer
    assert "separate summary and feature-CV workbooks." in result.answer


@pytest.mark.parametrize("message", ["how many proteins are in the plot?", "how was this plot made for my thesis?",
                                     "how do I save my session?", "write me a script to export this plot"])
def test_other_questions_are_not_download_questions(message):
    request = AssistantRequest(message=message, model="synthetic", view=_view())
    assert assistant._basic_download_response(request, 1) is None


def test_replicate_quality_on_a_qc_page_is_the_qc_review(assistant_rig):
    plan = json.dumps({"task": "current_plot"})
    assistant_rig.run([plan, "{}", "{}", "{}"], message="are my replicates ok?", view=_view())
    assert assistant_rig.plans[0].task == "qc_review"
    assert not any("sources" in call for call in assistant_rig.calls)


def _volcano_view(level="protein"):
    def select(id_, label, value, options):
        return {"id": id_, "label": label, "type": "select", "value": value, "options": options,
                "optionCount": len(options), "permission": "view"}
    return _view("stats.volcano", controls=[select("lvl", "Dataset level", level, ["protein", "phospho"]),
                                             select("c1", "Condition 1", "WT", ["WT", "GRK"]),
                                             select("c2", "Condition 2", "WT", ["WT", "GRK"])])


def test_the_phosphosite_level_is_set_before_the_comparison(assistant_rig):
    message = "volcano of the phosphosites, GRK vs WT"
    first = assistant_rig.run([], message=message, view=_volcano_view(), receipts=["Opened Statistical Analysis / Volcano Plot."])
    assert [(a.target, a.value) for a in first.actions] == [("lvl", "phospho")] and not first.finished
    second = assistant_rig.run([], message=message, view=_volcano_view("phospho"),
                               receipts=["Opened Statistical Analysis / Volcano Plot.", "Set Dataset level to Phospho"])
    assert [(a.target, a.value) for a in second.actions] == [("c1", "GRK")]
    assert assistant_rig.calls == []


def test_a_question_about_sites_does_not_switch_the_level(assistant_rig):
    request = AssistantRequest(message="how many sites are significant?", model="synthetic", view=_volcano_view())
    assert assistant._basic_level_plan(request) is None
