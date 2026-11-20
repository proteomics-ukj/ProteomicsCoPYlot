import { createContext, useContext, useState, type ComponentType } from "react";
import { readPersistedUiState } from "./persistedUiState";

export type Workflow = "protein" | "phospho" | "phosprot" | "peptide" | "precursor";
const workflows: Workflow[] = ["protein", "phospho", "phosprot", "peptide", "precursor"];

const WorkflowContext = createContext<{ kind: Workflow; select: (kind: Workflow) => void } | null>(null);

export function useWorkflow<K extends Workflow>(): [K, (kind: K) => void] {
  const scope = useContext(WorkflowContext);
  if (!scope) throw new Error("Missing workflow settings scope");
  return [scope.kind as K, scope.select];
}

// Remount the controls and cancel their pending requests together on a workflow switch.
export function withWorkflowSettings<P extends object>(
  key: string, Page: ComponentType<P>, available?: (props: P) => { value: string }[],
) {
  return function WorkflowPage(props: P) {
    const [preferred, select] = useState<Workflow>(() => {
      const saved = readPersistedUiState<{ kind?: Workflow }>(key, {}).kind;
      return saved && workflows.includes(saved) ? saved : "protein";
    });
    const choices = available?.(props);
    const kind = choices?.length && !choices.some(option => option.value === preferred)
      ? choices[0].value as Workflow : preferred;
    return <WorkflowContext.Provider value={{ kind, select }}><Page key={kind} {...props} /></WorkflowContext.Provider>;
  };
}
