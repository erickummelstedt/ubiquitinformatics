// Placeholder for API calls to backend (FastAPI)
export const fetchInitialReactions = async () => {
  const response = await fetch('/initialReactions.json');
  return response.json();
};
