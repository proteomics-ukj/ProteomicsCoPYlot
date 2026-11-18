import { useEffect } from "react";

export function useVisibleTab<T extends string>(
  activeTab: T, tabs: readonly { key: T }[], fallback: T, setActiveTab: (tab: T) => void,
) {
  useEffect(() => {
    if (!tabs.some((tab) => tab.key === activeTab)) {
      setActiveTab(tabs[0]?.key ?? fallback);
    }
  }, [activeTab, tabs, fallback, setActiveTab]);
}
