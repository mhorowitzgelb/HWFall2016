import sys




def find_aug_path(g):
    q = [(None,'r')]
    visited = {}

    while len(q) != 0:
        arc, node = q.pop(0)
        visited[node] = arc
        if node == 's':
            return visited
        for head in g[node]:
            if not head in visited:
                cap, flow = g[node][head]
                if flow < cap:
                    q.append(((node,head),head))
        for tail in g:
            if node in tail and not tail in visited:
                cap, flow = g[tail][node]
                if flow > 0:
                    q.append(((tail,node),tail))
    return visited


def find_flow(g):

    while 1:
        visited = find_aug_path(g)
        cur = 's'
        if not 's' in visited:
            break

        aug_path = []
        width = 999999
        while 1:
            if visited[cur] == None:
                break
            (head, tail) = visited[cur]
            (cap,flow) = g[head][tail]
            if cur == head:
                cur = tail
                aug_path.insert(0,(-1,(head,tail)))
                width = min(width,flow)
            else:
                cur = head
                aug_path.insert(0,(1,(head,tail)))
                width = min(width, cap-flow)
        print("Corrected augmented path")
        for (mult, (head,tail)) in aug_path:
            print(head+tail, end=" ")
            (cap,flow) = g[head][tail]
            flow += mult * width
            g[head][tail] = (cap,flow)
        print('width:', width)
    print("Final flow:")
    for head in g:
        for tail in g[head]:
            (cap,flow) = g[head][tail]
            print(head+tail, 'u:', cap, "x:", flow )

g = {
    'r' : {'p' : (6,0), 'a' : (9,0), 'q' : (4,0)},
    'p' : {'q' : (2,0), 'b' : (3,0)},
    'q' : {'p' : (1,0), 'b' : (2,0), 'd' : (6,0)},
    'a' : {'c': (8,0), 'd' : (1,0)},
    'b' : {'a' : (1,0), 's' : (8,0)},
    'c' : {'b' : (2,0), 's' : (4,0) , 'q' : (1,0)},
    'd' : {'c' : (1,0), 's' : (6,0)}
}

find_flow(g)




