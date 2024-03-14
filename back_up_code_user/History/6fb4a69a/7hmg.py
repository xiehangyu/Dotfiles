import pygame
from queue import LifoQueue
import copy
# Initialize Pygame
def simplemanualsimplification(size1,size2,initial_pos=0,final_pos=0,supportlength=0):
    pygame.init()
    lifo_queue=LifoQueue()
    # Set the dimensions of the window
    win_width = 0
    win_height = 0
    win = pygame.display.set_mode((win_width, win_height))


    # Define the brick size and color
    brick_size = 50
    brick_color = (0, 255, 0)  # Green

    # Generate the initial set of bricks in brickwall pattern
    bricks = [(2 * i + (2 * ((j + 1) // 2 - j / 2)), j) for j in range(size2) for i in range(size1)]
    bottomstartingpointx=(2*((size2)//2-(size2-1)/2))
    bottomstartingpointy=size2
    circles=[(i,0)for i in range(final_pos,final_pos+supportlength)]+[(bottomstartingpointx+i,bottomstartingpointy) for i in range(initial_pos,initial_pos+supportlength)]
    print(circles)
    # Game loop
    run = True
    while run:
        pygame.time.delay(100)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_LEFT:
                    print(event.key)
                    print("haha")
                    if lifo_queue.empty()==False:
                        bricks=lifo_queue.get()
                        print(lifo_queue.qsize())
                        print(bricks)
            elif event.type == pygame.MOUSEBUTTONDOWN:
                # Get the position of the click
                mouse_pos = pygame.mouse.get_pos()

                # Calculate the coordinates of the brick that was clicked
                clicked_brick = (mouse_pos[0] // brick_size, mouse_pos[1] // brick_size)

                # Remove the brick from the list
                if clicked_brick in bricks:
                    lifo_queue.put(copy.deepcopy(bricks))
                    bricks.remove(clicked_brick)

        # Clear the window
        win.fill((255, 255, 255))

        # Draw the bricks
        for brick in bricks:
            x, y = brick
            pygame.draw.rect(win, brick_color, (x * brick_size, y * brick_size, brick_size, brick_size))
            pygame.draw.line(win, (0, 0, 0), (x * brick_size, y * brick_size),
                             ((x + 1) * brick_size, (y + 1) * brick_size), 3)
            pygame.draw.line(win, (0, 0, 0), ((x + 1) * brick_size, y * brick_size),
                             (x * brick_size, (y + 1) * brick_size), 3)
        for circle in circles:
            x,y=circle
            pygame.draw.circle(win, (0, 0, 0), (x * brick_size, y * brick_size), 10)

        # Update the window
        pygame.display.update()

    pygame.quit()



if __name__ == '__main__':
    simplemanualsimplification(6,5,4,7,1)