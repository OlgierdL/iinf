'''Functions to enumerate all perfect and maximum matchings in bipartited graph.
Implemented following the algorithms in the paper "Algorithms for Enumerating All Perfect, Maximum and Maximal Matchings in Bipartite Graphs" by Takeaki Uno,
using numpy and networkx modules of python.

Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
'''
#--------Import modules-------------------------
import networkx as nx
from networkx import bipartite
import numpy
from scipy import sparse
import sys

def formDirected(g,match):
    '''Form directed graph D from G and matching M.

    <g>: undirected bipartite graph. Nodes are separated by their
         'bipartite' attribute.
    <match>: list of edges forming a matching of <g>.

    Return <d>: directed graph, with edges in <match> pointing from set-0
                (bipartite attribute ==0) to set-1 (bipartite attrbiute==1),
                and the other edges in <g> but not in <matching> pointing
                from set-1 to set-0.
    '''
    d=nx.DiGraph()
    for ee in g.edges():
        if ee in match or (ee[1],ee[0]) in match:
            if g.nodes[ee[0]]['bipartite']==0:
                d.add_edge(ee[0],ee[1])
            else:
                d.add_edge(ee[1],ee[0])
        else:
            if g.nodes[ee[0]]['bipartite']==0:
                d.add_edge(ee[1],ee[0])
            else:
                d.add_edge(ee[0],ee[1])
    return d

def enumMaximumMatching(g, m=sys.maxsize):
    '''Find all maximum matchings in an undirected bipartite graph.

    <g>: undirected bipartite graph. Nodes are separated by their
         'bipartite' attribute.
    <m>: maximum number of matchings to find.

    Return <all_matches>: list, each is a list of edges forming a maximum
                          matching of <g>.

    Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
    Update time: 2017-05-21 20:04:51.
    '''
    all_matches=[]
    if m < 1:
        return all_matches
    #----------------Find one matching M----------------
    #nodes = g.nodes
    top = [n for n in g.nodes if g.nodes[n]['bipartite']==0]
    match=bipartite.hopcroft_karp_matching(g, top_nodes=top)
    #---------------Re-orient match arcs---------------
    match2=[]
    for kk,vv in match.items():
        if g.nodes[kk]['bipartite']==0:
            match2.append((kk,vv))
    match=match2
    all_matches.append(match)
    if len(all_matches) == m:
        return all_matches
    #-----------------Enter recursion-----------------
    all_matches=enumMaximumMatchingIter(g,match,all_matches,m,None)
    return all_matches

def enumMaximumMatchingIter(g,match,all_matches,m,add_e=None):
    '''Recurively search maximum matchings.

    <g>: undirected bipartite graph. Nodes are separated by their
         'bipartite' attribute.
    <match>: list of edges forming one maximum matching of <g>.
    <all_matches>: list, each is a list of edges forming a maximum
                   matching of <g>. Newly found matchings will be appended
                   into this list.
    <add_e>: tuple, the edge used to form subproblems. If not None,
             will be added to each newly found matchings.

    Return <all_matches>: updated list of all maximum matchings.

    Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
    Update time: 2017-05-21 20:09:06.
    '''
    #---------------Form directed graph D---------------
    d=formDirected(g,match)
    #-----------------Find cycles in D-----------------
    try:
        cycle=nx.find_cycle(d)
    except nx.exception.NetworkXNoCycle:
        cycle=[]
    if len(cycle)==0:
        #---------If no cycle, find a feasible path---------
        all_uncovered=set(g.nodes).difference(set([ii[0] for ii in match]))
        all_uncovered=all_uncovered.difference(set([ii[1] for ii in match]))
        all_uncovered=list(all_uncovered)
        #--------------If no path, terminiate--------------
        if len(all_uncovered)==0:
            return all_matches
        #----------Find a length 2 feasible path----------
        idx=0
        uncovered=all_uncovered[idx]
        while True:
            if uncovered not in nx.isolates(g):
                paths=nx.single_source_shortest_path(d,uncovered,cutoff=2)
                len2paths=[vv for kk,vv in paths.items() if len(vv)==3]
                if len(len2paths)>0:
                    reversed=False
                    break
                #----------------Try reversed path----------------
                paths_rev=nx.single_source_shortest_path(d.reverse(),uncovered,cutoff=2)
                len2paths=[vv for kk,vv in paths_rev.items() if len(vv)==3]
                if len(len2paths)>0:
                    reversed=True
                    break
            idx+=1
            if idx>len(all_uncovered)-1:
                return all_matches
            uncovered=all_uncovered[idx]
        #-------------Create a new matching M'-------------
        len2path=len2paths[0]
        if reversed:
            len2path=len2path[::-1]
        len2path=list(zip(len2path[:-1],len2path[1:]))
        new_match=[]
        for ee in d.edges():
            if ee in len2path:
                if g.nodes[ee[1]]['bipartite']==0:
                    new_match.append((ee[1],ee[0]))
            else:
                if g.nodes[ee[0]]['bipartite']==0:
                    new_match.append(ee)
        if add_e is not None:
            for ii in add_e:
                new_match.append(ii)
        all_matches.append(new_match)
        if len(all_matches) == m:
            return all_matches
        #---------------------Select e---------------------
        e=set(len2path).difference(set(match))
        e=list(e)[0]
        #-----------------Form subproblems-----------------
        g_plus=g.copy()
        g_minus=g.copy()
        g_plus.remove_node(e[0])
        g_plus.remove_node(e[1])
        g_minus.remove_edge(e[0],e[1])
        add_e_new=[e,]
        if add_e is not None:
            add_e_new.extend(add_e)
        all_matches=enumMaximumMatchingIter(g_minus,match,all_matches,m,add_e)
        if len(all_matches) == m:
            return all_matches
        all_matches=enumMaximumMatchingIter(g_plus,new_match,all_matches,m,add_e_new)
    else:
        #-------------Create a new matching M'-------------
        new_match=[]
        for ee in d.edges():
            if ee in cycle:
                if g.nodes[ee[1]]['bipartite']==0:
                    new_match.append((ee[1],ee[0]))
            else:
                if g.nodes[ee[0]]['bipartite']==0:
                    new_match.append(ee)
        if add_e is not None:
            for ii in add_e:
                new_match.append(ii)
        all_matches.append(new_match)
        if len(all_matches) == m:
            return all_matches
        #-----------------Choose an edge E-----------------
        e=set(match).intersection(set(cycle))
        e=list(e)[0]
        #-----------------Form subproblems-----------------
        g_plus=g.copy()
        g_minus=g.copy()
        g_plus.remove_node(e[0])
        g_plus.remove_node(e[1])
        g_minus.remove_edge(e[0],e[1])
        add_e_new=[e,]
        if add_e is not None:
            add_e_new.extend(add_e)
        all_matches=enumMaximumMatchingIter(g_minus,new_match,all_matches,m,add_e)
        if len(all_matches) == m:
            return all_matches
        all_matches=enumMaximumMatchingIter(g_plus,match,all_matches,m,add_e_new)
    return all_matches

def compute_all_matches(edges, tcount, mcount):
    g=nx.Graph()
    for i in range(tcount):
        g.add_node(i, bipartite=0)
    for i in range(mcount):
        g.add_node(i+tcount,bipartite=1)
    g.add_edges_from(edges)
    return enumMaximumMatching(g)

def example():
    g=nx.Graph()
    n = 4
    edges=[
            [0, 4],[0, 5],
            [1, 6],[1, 7],
            [2, 5],[2, 7],
            [3, 4],[3, 6]
            ]
    for i in range(n):
        g.add_node(i,   bipartite=0) 
        g.add_node(i+n, bipartite=1)
    g.add_edges_from(edges)
    all_matches=enumMaximumMatching(g)
    print(all_matches)

#-------------Main---------------------------------
if __name__=='__main__':
    example()