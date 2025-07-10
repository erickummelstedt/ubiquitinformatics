import { simulateClicksFromJson } from '../src/components/ScaffoldJsonWrapper';

const DEFAULT_NODES = [
  { x: 300, y: 300, clicks: 0 },
  { x: 250, y: 250, clicks: 0 }, { x: 350, y: 250, clicks: 0 },
  { x: 200, y: 200, clicks: 0 }, { x: 300, y: 200, clicks: 0 }, { x: 400, y: 200, clicks: 0 },
  { x: 150, y: 150, clicks: 0 }, { x: 250, y: 150, clicks: 0 }, { x: 350, y: 150, clicks: 0 }, { x: 450, y: 150, clicks: 0 },
  { x: 100, y: 100, clicks: 0 }, { x: 200, y: 100, clicks: 0 }, { x: 300, y: 100, clicks: 0 }, { x: 400, y: 100, clicks: 0 }, { x: 500, y: 100, clicks: 0 }
];

const DEFAULT_EDGES = [
  [0, 1], [0, 2],
  [1, 3], [1, 4],
  [2, 4], [2, 5],
  [3, 6], [3, 7],
  [4, 7], [4, 8],
  [5, 8], [5, 9],
  [6, 10], [6, 11],
  [7, 11], [7, 12],
  [8, 12], [8, 13],
  [9, 13], [9, 14]
];

const TEST_JSON = {
  chain_number: 1,
  branching_sites: [
    {
      site_name: 'K48',
      children: {
        chain_number: 2,
        branching_sites: [
          {
            site_name: 'K63',
            children: {
              chain_number: 3,
              branching_sites: []
            }
          }
        ]
      }
    },
    {
      site_name: 'K63',
      children: {
        chain_number: 4,
        branching_sites: []
      }
    }
  ]
};

test('simulateClicksFromJson generates correct nodes and arrows', () => {
  const { nodes, arrows } = simulateClicksFromJson(TEST_JSON, DEFAULT_NODES.map(n => ({ ...n })), DEFAULT_EDGES);

  // Check node clicks
  expect(nodes[0].clicks).toBe(1); // Node 1
  expect(nodes[2].clicks).toBe(1); // Node 2
  expect(nodes[3].clicks).toBe(1); // Node 3
  expect(nodes[1].clicks).toBe(1); // Node 4

  // Check arrows
  expect(arrows).toEqual([
    { from: 0, to: 2, color: 'red', highlighted: true },
    { from: 2, to: 3, color: 'blue', highlighted: true },
    { from: 0, to: 1, color: 'blue', highlighted: true }
  ]);
});
