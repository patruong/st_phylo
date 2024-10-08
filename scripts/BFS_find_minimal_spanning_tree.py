from collections import deque

def calculate_distance(v1, v2):
    return sum(abs(a - b) for a, b in zip(v1, v2))

def create_graph(vectors, max_distance=1.0):
    graph = {}
    for v in vectors:
        graph[tuple(v)] = set()
    
    for v1 in vectors:
        for v2 in vectors:
            if v1 != v2:
                distance = calculate_distance(v1, v2)
                if distance <= max_distance:  # Adjust max_distance to allow larger jumps
                    graph[tuple(v1)].add(tuple(v2))
                    graph[tuple(v2)].add(tuple(v1))
    
    return graph

def bfs_shortest_path(graph, start, target, heuristic_threshold=0.5):
    queue = deque([[start]])
    visited = set([start])
    
    while queue:
        path = queue.popleft()
        node = path[-1]
        
        if node == target:
            return path
        
        # Sort neighbors by their distance from the current node
        neighbors = sorted(graph[node], key=lambda neighbor: calculate_distance(node, neighbor))
        
        for neighbor in neighbors:
            distance = calculate_distance(node, neighbor)
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
                # Only prioritize smaller changes (heuristic threshold) but allow larger jumps if needed
                if distance <= heuristic_threshold:
                    break  # Prioritize the closest neighbor by breaking the loop early
    
    return None  # No path found

# Example usage
vectors = [
    [1.0, 1.0, 1.0, 1.0],
    [1.0, 1.5, 1.0, 2.0],
    [1.0, 1.0, 2.0, 1.5],
    [1.0, 1.0, 1.0, 1.5],
    [1.0, 1.5, 1.0, 1.5],
    [1.0, 1.5, 1.0, 2.0],
    [1.0, 1.5, 1.0, 3.0]
]

# Adjust the max_distance to include larger jumps if needed
graph = create_graph(vectors, max_distance=2.0)

start = tuple([1.0, 1.0, 1.0, 1.0])
target1 = tuple([1.0, 1.5, 1.0, 3.0])
target2 = tuple([1.0, 1.0, 2.0, 1.5])

# Run BFS with a heuristic threshold of 0.5 for small changes
path1 = bfs_shortest_path(graph, start, target1, heuristic_threshold=0.5)
path2 = bfs_shortest_path(graph, start, target2, heuristic_threshold=0.5)

print("Shortest path to [1.0, 1.5, 1.0, 2.0]:")
for step in path1:
    print(step)

print("\nShortest path to [1.0, 1.0, 2.0, 1.5]:")
for step in path2:
    print(step)

print("\nPath distances:")
for i, path in enumerate([path1, path2], 1):
    if path:
        total_distance = sum(calculate_distance(path[j], path[j+1]) for j in range(len(path)-1))
        print(f"Path {i}: {total_distance}")
    else:
        print(f"Path {i}: No path found")


