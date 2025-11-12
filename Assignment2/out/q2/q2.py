import networkx as nx

G = nx.Graph()
G.add_edges_from([("A", "x"), ("B", "x"), ("C", "x"),
                 ("D", "x"), ("E", "x"), ("D", "E")])

betweenness = nx.betweenness_centrality(G)

print(betweenness)
