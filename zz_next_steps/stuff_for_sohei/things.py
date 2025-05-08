import pygame

# Pygame setup
pygame.init()
screen = pygame.display.set_mode((800, 600))
pygame.display.set_caption("Ubiquitin-like Graph Builder with Base Structure")
clock = pygame.time.Clock()

# Colors
WHITE = (255, 255, 255)
GRAY = (180, 180, 180)
LIGHT_GRAY = (220, 220, 220)
BLACK = (0, 0, 0)

radius = 20

# Base structure: 5 complete levels with cleaner spacing
nodes = [
    # Level 0 (root)
    (400, 500),  # 0

    # Level 1
    (300, 400),  # 1
    (500, 400),  # 2

    # Level 2
    (200, 300),  # 3
    (400, 300),  # 4
    (600, 300),  # 5

    # Level 3
    (100, 200),  # 7
    (300, 200),  # 8
    (500, 200),  # 9
    (700, 200),  # 10

    # Level 4
    (0, 100),  # 13
    (200, 100),  # 14
    (400, 100),  # 16
    (600, 100),  # 18
    (800, 100)  # 19
]

edges = [
    # Level 0 to Level 1
    (0, 1), (0, 2),

    # Level 1 to Level 2
    (1, 3), (1, 4),
    (2, 4), (2, 5),

    # Level 2 to Level 3
    (3, 6), (3, 7),
    (4, 7), (4, 8),
    (5, 8), (5, 9),

    # Level 3 to Level 4
    (6, 10), (6, 11),
    (7, 11), (7, 12),
    (8, 12), (8, 13),
    (9, 13), (9, 14)
]

# For new additions
last_node_idx = None

def draw():
    screen.fill(WHITE)

    # Draw edges
    for i1, i2 in edges:
        pygame.draw.line(screen, GRAY, nodes[i1], nodes[i2], 2)

    # Draw nodes
    for idx, (x, y) in enumerate(nodes):
        pygame.draw.circle(screen, LIGHT_GRAY, (x, y), radius)
        pygame.draw.circle(screen, BLACK, (x, y), radius, 2)
        font = pygame.font.SysFont(None, 24)
        img = font.render(str(idx), True, BLACK)
        screen.blit(img, (x - 10, y - 10))

    pygame.display.flip()

running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

        elif event.type == pygame.MOUSEBUTTONDOWN:
            x, y = pygame.mouse.get_pos()
            nodes.append((x, y))
            current_idx = len(nodes) - 1

            if last_node_idx is not None:
                lx, ly = nodes[last_node_idx]
                x, y = lx, ly - 80  # Go upward
            else:
                x, y = pygame.mouse.get_pos()
                
            last_node_idx = current_idx

    draw()
    clock.tick(60)

pygame.quit()