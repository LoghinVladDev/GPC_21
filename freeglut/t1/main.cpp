#include <GL/glut.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>

unsigned char prevKey;

void Display1() {
    glColor3f(0.2,0.15,0.88);
    glBegin(GL_LINES);
    glVertex2i(1,1);
    glVertex2i(-1,-1);
    glEnd();

    glColor3f(1,0.1,0.1);
    glBegin(GL_LINES);
    glVertex2i(-1,1);
    glVertex2i(1,-1);
    glEnd();

    glBegin(GL_LINES);
    glVertex2d(-0.5,0);
    glVertex2d(0.5,0);
    glEnd();
}

void Display2() {
    glColor3f(1,0.1,0.1);
    glBegin(GL_LINES);
    glVertex2f(1.0,1.0);
    glVertex2f(0.9,0.9);
    glVertex2f(0.8,0.8);
    glVertex2f(0.7,0.7);
    glVertex2f(0.6,0.6);
    glVertex2f(0.5,0.5);
    glVertex2f(-0.5,-0.5);
    glVertex2f(-1.0,-1.0);
    glEnd();
}

void Display3() {
    // trasare puncte GL_POINTS : deseneaza n puncte
    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_POINTS);
    glVertex2f(0.5f, 0.5f);
    glVertex2f(0.5f, -0.5f);
    glVertex2f(-0.5f, 0.5f);
    glVertex2f(-0.5f, -0.5f);
    glEnd();
}

void Display4() {
    glColor3f(1,0.1,0.1); // rosu
    // trasare linie poligonala GL_LINE_STRIP : (v0,v1), (v1,v2), (v_{n-2},v_{n-1})
    glBegin(GL_LINE_STRIP);
    // de completat ...
    glVertex2f(1.0f, 1.0f);
    glVertex2f(1.0f, 0.85f);
    glVertex2f(0.6f, 0.65f);
    glVertex2f(0.6f, 0.5f);
    glEnd();
}

void Display5() {
    glColor3f(1,0.1,0.1); // rosu
    // trasare linie poligonala inchisa GL_LINE_LOOP : (v0,v1), (v1,v2), (v_{n-1},v0)
    glBegin(GL_LINE_LOOP);
    // de completat ...
    glVertex2f(1.0f, 1.0f);
    glVertex2f(1.0f, 0.85f);
    glVertex2f(0.6f, 0.65f);
    glVertex2f(0.6f, 0.5f);
    glEnd();
}

void Display6() {
    glPolygonMode(GL_FRONT, GL_FILL);
    glColor3f(1,0.1,0.1); // rosu
    // trasare triunghiuri GL_TRIANGLES : (v0,v1,v2), (v3,v4,v5), ...
    glBegin(GL_TRIANGLES);
    // de completat ...
    glVertex2f(-1.0f, -1.0f);
    glVertex2f(-0.8f, -0.8f);
    glVertex2f(-1.0f, -0.8f);

    glVertex2f(0.8f, 0.8f);
    glVertex2f(1.0f, 0.8f);
    glVertex2f(1.0f, 1.0f);
    glEnd();
}

void Display7() {
    // trasare patrulatere GL_QUADS : (v0,v1,v2,v3), (v4,v5,v6,v7), ...
    glColor3f(1,0.1,0.1); // rosu
    glPolygonMode(GL_FRONT, GL_FILL);
    glBegin(GL_QUADS);
    // de completat ...
    glVertex2f(1.0f, 1.0f);
    glVertex2f(0.3f, 0.7f);
    glVertex2f(0.5f, 0.5f);
    glVertex2f(1.0f, 0.5f);
    glEnd();
}

#include <iostream>
void Display8() {
    struct Point { float x; float y; Point(float x, float y) : x(x), y(y) {} };
    // trasare poligon convex GL_QUADS : (v0,v1,v2, ..., v_{n-1})
//    glColor3f(1,0.1,0.1); // rosu
//    glColor3f(0.2,0.15,0.88);
    Point centre = {0.0f, 0.0f};
    constexpr auto pi = static_cast <float> (M_PI);

    auto getHexPoint = [&centre] (unsigned int pointNumber, float dFromCentre) -> Point {
        pointNumber = pointNumber % 6;
        auto degrees = static_cast<float>(pointNumber * 60);

        return {
            centre.x + cosf ( degrees * pi / 180.0f ) * dFromCentre,
            centre.y + sinf ( degrees * pi / 180.0f ) * dFromCentre
        };
    };

    auto drawHex = [& getHexPoint] (float r, float g, float b, float radius) -> void {
        glColor3f(r,g,b);
        glBegin(GL_POLYGON);

        for ( int i = 0; i < 6; i++ ) {
            auto point = getHexPoint ( i, radius );
            glVertex2f( point.x, point.y );
        }

        glEnd();
    };

    glPolygonMode(GL_FRONT, GL_FILL);

    // de completat ...
    drawHex(0.2f, 0.15f, 0.88f, 0.75f);
    drawHex(1.0f, 0.1f, 0.1f, 0.5f);
    drawHex(1.0f, 1.0f, 1.0f, 0.48f);
}

void Init(void) {
    glClearColor(1.0,1.0,1.0,1.0);

    glLineWidth(3);
    glPointSize(4);

    // functia void glPolygonMode (GLenum face, GLenum mode)
    // controleaza modul de desenare al unui poligon
    // mode : GL_POINT (numai vf. primitivei) GL_LINE (numai
    //        muchiile) GL_FILL (poligonul plin)
    // face : tipul primitivei geometrice dpdv. al orientarii
    //        GL_FRONT - primitive orientate direct
    //        GL_BACK  - primitive orientate invers
    //        GL_FRONT_AND_BACK  - ambele tipuri
    glPolygonMode(GL_FRONT, GL_LINE);
}

void Display(void) {
    printf("Call Display\n");

    glClear(GL_COLOR_BUFFER_BIT);

    switch(prevKey) {
        case '1':
            Display1();
            break;
        case '2':
            Display2();
            break;
        case '3':
            Display3();
            break;
        case '4':
            Display4();
            break;
        case '5':
            Display5();
            break;
        case '6':
            Display6();
            break;
        case '7':
            Display7();
            break;
        case '8':
            Display8();
            break;
        default:
            break;
    }

    glFlush();
}

void Reshape(int w, int h) {
    printf("Call Reshape : latime = %d, inaltime = %d\n", w, h);
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}

void KeyboardFunc(unsigned char key, int x, int y) {
    printf("Ati tastat <%c>. Mouse-ul este in pozitia %d, %d.\n",key, x, y);
    prevKey = key;
    if (key == 27)
        exit(0);
    glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) {
    printf("Call MouseFunc : ati %s butonul %s in pozitia %d %d\n",
        (state == GLUT_DOWN) ? "apasat" : "eliberat",
        (button == GLUT_LEFT_BUTTON) ? "stang" : ((button == GLUT_RIGHT_BUTTON) ? "drept": "mijlociu"),
        x, y);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(300, 300);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutCreateWindow (argv[0]);

    Init();

    glutReshapeFunc(Reshape);
    glutKeyboardFunc(KeyboardFunc);
    glutMouseFunc(MouseFunc);
    glutDisplayFunc(Display);
    glutMainLoop();

    return 0;
}