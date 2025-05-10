 
import turtle

turtle.hideturtle()
turtle.tracer(0)
turtle.penup()
turtle.setpos(-100, -150)
turtle.pendown()

axiom, tempAx, logic, count = 'FX', '', {'X': 'X+YF+', 'Y': '-FX-Y'}, 200

for i in range(count):
    for j in axiom:
        tempAx += logic[j] if j in logic else j
    axiom, tempAx = tempAx, ''

for k in axiom:
    if k == 'F':
        turtle.forward(2.5)
    elif k == '+':
        turtle.right(90)
    elif k == '−':
        turtle.left(90)

turtle.update()
turtle.mainloop()