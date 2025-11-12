import networkx as nx
import numpy as np


def main():
    G = nx.barabasi_albert_graph(n=30, m=2, seed=42)

    A = nx.to_numpy_array(G, dtype=int, nodelist=sorted(G.nodes()))
    np.savetxt("adjacency_matrix.csv", A, fmt="%d", delimiter=",")

    nx.write_graphml(G, "network.graphml")


if __name__ == "__main__":
    main()
