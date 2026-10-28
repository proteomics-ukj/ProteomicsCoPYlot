"""Public assistant errors: return approved text, never arbitrary exception details."""
import re


DEFAULT_ERROR = "The assistant could not complete this step. Check your input and AI connection settings, then retry."

_MESSAGES = (
    "The assistant request was cancelled.",
    "The assistant time limit was reached.",
    "Provider operation exceeded its time limit.",
    "Provider response exceeds the configured size limit.",
    "The analysis changed while the assistant was working. Please retry.",
    "This request exceeds the remaining action budget. Open the canvas and ask for a smaller batch.",
    "The analysis canvas is not available in this view.",
    "The analysis canvas is unavailable in this view.",
    "A feature-focused canvas needs one explicit feature scope and no other actions.",
    "Expanding a protein into all sites requires an explicit all-sites request.",
    "The requested data level must not be substituted.",
    "Enter the address of the local AI server, for example http://127.0.0.1:8080/v1.",
    "The local AI server address must be an http or https URL.",
    "The local AI server address must not contain credentials; use the API key field.",
    "The local AI server address must not contain a query string or fragment.",
    "The local AI server address has an invalid port.",
    "The local AI server redirected the request. Configure the final address directly.",
    "Could not reach the local AI server. Check the address and connection.",
    "The local AI server refused the request; check its API key.",
    "No model is selected for the local AI server.",
    "Generated AI context contains a protected path or credential-like value and was not sent.",
    "The assistant never changes how data are filtered, normalized, imputed or transformed, "
    "or the localization cutoff. Use the analysis app to change processing settings manually. "
    "No data-processing settings were changed.",
) + tuple(
    f"{label} API key is required. Set {env} on the server, or enter a key in the chat settings."
    for label, env in (("OpenAI", "OPENAI_API_KEY"), ("Anthropic", "ANTHROPIC_API_KEY"),
                       ("Google Gemini", "GEMINI_API_KEY"))
) + tuple(
    message
    for provider in ("OpenAI", "Anthropic", "Gemini", "Ollama", "The local AI server")
    for message in (
        f"Could not reach {provider}. Check the endpoint and connection.",
        f"{provider} redirect blocked. Configure the final endpoint directly.",
        *(f"{provider} request failed (HTTP {status})." for status in range(400, 600)),
    )
)


def public_assistant_error(error: Exception) -> str:
    """Match complete messages and return our own text, not exception-derived strings.

    Unknown library errors may contain file paths, input snippets or credentials.
    Prefix/substring matches must not allow such suffixes through this boundary.
    """
    if isinstance(error, ValueError) and len(error.args) == 1 and isinstance(error.args[0], str):
        for message in _MESSAGES:
            if error.args[0] == message:
                return message
        if re.fullmatch(
            r"(?:OpenAI|Anthropic|Gemini|Ollama) did not answer within [0-9]+ seconds; "
            r"this step was stopped\. Earlier completed actions are retained\.", error.args[0]
        ):
            return "The AI provider did not answer in time. Earlier completed actions are retained. Please retry."
        if re.fullmatch(
            r"The local AI server did not answer within [0-9]+ seconds, so this step was stopped\. "
            r"Earlier completed actions are retained\. Try again, ask for one thing at a time, "
            r"or raise COPYLOT_LOCAL_AI_TIMEOUT if the server is slow\.", error.args[0]
        ):
            return ("The local AI server did not answer in time. Earlier completed actions are retained. "
                    "Try again, ask for one thing at a time, or raise COPYLOT_LOCAL_AI_TIMEOUT if the server is slow.")
    return DEFAULT_ERROR
