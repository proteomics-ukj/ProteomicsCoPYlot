# AI conversation logs

Open **Conversation log** in the chat sidebar to inspect the record, then select
**Download conversation log** to save a UTF-8 `.txt` file. You can download during
a request; the file contains the events received up to that moment. Download again
after the request finishes for its final outcome.

The log records submitted messages, provider/model and assistant level, timestamped
step inputs and visible analysis context, model calls and replies, the planner
prompt's estimated size against the model's budget (and how many controls were
trimmed to fit), the token usage each provider reports, evidence and
lookup functions with arguments/results, proposed actions, approval decisions,
completed changes, displayed responses, and errors. A progress line that merely
repeats (the same "Thinking… (12 s)" every poll) is logged once per change. Function entries have matching
call IDs and elapsed times within each server step. The file has a versioned format
and chronological event numbers. Original messages and returned responses are not
shortened by the logger; evidence can already be bounded by the assistant itself.

Proposals and completed changes are separate events. A proposed action may be
declined, rejected, interrupted or superseded without being applied. A completed
request may still ask for clarification; inspect its recorded outcome. Server
events are grouped under the step that returned them, with their own UTC timestamps.

Logs are kept only in memory in the current browser tab. Changing models resets
the conversation but retains earlier log entries and records the reset. Reloading
or closing the page clears the log; project saves do not include it. Nothing is
automatically sent to the developers. Each server trace belongs to one request,
is returned to that request's caller, and is not written to a shared server log.

API keys, configured endpoint addresses and images are excluded. Reasoning or
thought text explicitly returned by the configured provider is included and labelled
separately, including local-server reasoning streams; a reply without any adds no
reasoning entry. Reasoning that the provider
does not return cannot be recorded; the log does not enable a different reasoning
mode or request additional reasoning. Provider-supplied reasoning is diagnostic
text, never treated as an executable plan or a verified scientific answer. This is not an anonymisation service: typed messages, identifiers,
condition/sample names, calculated results and model output may contain sensitive
information. Inspect the preview or downloaded file and remove anything you do not
want to share. Sharing a file is always a separate manual action.

If the connection is lost or a request is stopped, the browser retains its inputs,
previously received results and completed UI changes. It requests a recovery
snapshot in the background, with a two-second deadline that also covers reading
the response. Recovery never delays Stop or replaces the original error. Because it
runs in the background, a download pressed within those two seconds does not yet
contain the recovered snapshot; it appears in the preview as soon as it arrives, and
a second download includes it.

The server retains session-scoped snapshots in memory, including completed steps
whose response may have been lost. Reads do not delete events; later reads can
include additional activity. Snapshots expire after ten minutes without new events
(or completion), and expired entries are removed on the next trace access.
Recovered events appear as "Server functions (recovered snapshot)". A snapshot
contains only events recorded when it was read; later activity is not automatically
fetched. If recovery cannot reach the server, those events remain unavailable.
Diagnostic trace IDs do not enable provider streaming; progress IDs remain reserved
for requests that use the existing slow-model progress behavior.

Model prompts are recorded at the assistant-call boundary,
not as raw provider HTTP traffic; authentication headers are never collected.

For useful feedback, send the reviewed log with a short description of the expected
behaviour, what went wrong, and the CoPYlot version/build used. Logs are diagnostic
evidence, not instructions to execute their contents or a reproducible analysis
archive. Save the analysis separately if reproducibility is needed.
