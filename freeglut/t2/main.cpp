#include <GL/glut.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// dimensiunea ferestrei in pixeli
#define dim 300

unsigned char prevKey;

// concoida lui Nicomede (concoida dreptei)
// $x = a + b \cdot cos(t), y = a \cdot tg(t) + b \cdot sin(t)$. sau
// $x = a - b \cdot cos(t), y = a \cdot tg(t) - b \cdot sin(t)$. unde
// $t \in (-\pi / 2, \pi / 2)$
void Display1() {
    double xmax, ymax, xmin, ymin;
    double a = 1, b = 2;
    double pi = 4 * atan(1);
    double ratia = 0.05;

    // calculul valorilor maxime/minime ptr. x si y
    // aceste valori vor fi folosite ulterior la scalare
    xmax = a - b - 1;
    xmin = a + b + 1;
    ymax = ymin = 0;
    for (double t = - pi/2 + ratia; t < pi / 2; t += ratia) {
        double x1, y1, x2, y2;
        x1 = a + b * cos(t);
        xmax = (xmax < x1) ? x1 : xmax;
        xmin = (xmin > x1) ? x1 : xmin;

        x2 = a - b * cos(t);
        xmax = (xmax < x2) ? x2 : xmax;
        xmin = (xmin > x2) ? x2 : xmin;

        y1 = a * tan(t) + b * sin(t);
        ymax = (ymax < y1) ? y1 : ymax;
        ymin = (ymin > y1) ? y1 : ymin;

        y2 = a * tan(t) - b * sin(t);
        ymax = (ymax < y2) ? y2 : ymax;
        ymin = (ymin > y2) ? y2 : ymin;
    }

    xmax = (fabs(xmax) > fabs(xmin)) ? fabs(xmax) : fabs(xmin);
    ymax = (fabs(ymax) > fabs(ymin)) ? fabs(ymax) : fabs(ymin);

    // afisarea punctelor propriu-zise precedata de scalare
    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (double t = - pi/2 + ratia; t < pi / 2; t += ratia) {
        double x1, y1, x2, y2;
        x1 = (a + b * cos(t)) / xmax;
        x2 = (a - b * cos(t)) / xmax;
        y1 = (a * tan(t) + b * sin(t)) / ymax;
        y2 = (a * tan(t) - b * sin(t)) / ymax;

        glVertex2f(x1,y1);
    }
    glEnd();

    glBegin(GL_LINE_STRIP);
    for (double t = - pi/2 + ratia; t < pi / 2; t += ratia) {
        double x1, y1, x2, y2;
        x1 = (a + b * cos(t)) / xmax;
        x2 = (a - b * cos(t)) / xmax;
        y1 = (a * tan(t) + b * sin(t)) / ymax;
        y2 = (a * tan(t) - b * sin(t)) / ymax;

        glVertex2f(x2,y2);
    }
    glEnd();
}

// graficul functiei
// $f(x) = \bar sin(x) \bar \cdot e^{-sin(x)}, x \in \langle 0, 8 \cdot \pi \rangle$,
void Display2() {
    double pi = 4 * atan(1);
    double xmax = 8 * pi;
    double ymax = exp(1.1);
    double ratia = 0.05;

    // afisarea punctelor propriu-zise precedata de scalare
    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_LINE_STRIP);
    for (double x = 0; x < xmax; x += ratia) {
        double x1, y1;
        x1 = x / xmax;
        y1 = (fabs(sin(x)) * exp(-sin(x))) / ymax;

        glVertex2f(x1,y1);
    }
    glEnd();
}

void Display3() {
    constexpr float xMin = 0.0f;
    constexpr float xMax = 100.0f;

    constexpr float step = 0.05f;

    constexpr auto nearestIntDistance = [] (float x) -> float {
        auto low = x - static_cast<float>(static_cast<int>(x));
        auto high = static_cast<float>(static_cast<int>(x) + 1) - x;

        return std::min(low, high);
    };

    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for ( float x = xMin; x <= xMax; x += step ) {
//        printf("MD : %lf\n", nearestIntDistance(x));
        float y = x == 0.0f ? 1.0f : ( nearestIntDistance(x) / x );

        glVertex2f( x / 30.0f, y );
    }

    glEnd();
}

void Display4() {
    constexpr float a = 0.3f, b = 0.2f;

    constexpr float step = 0.05f;
    constexpr auto pi = static_cast <float> (M_PI);
    constexpr auto minT = - pi;
    constexpr auto maxT = pi;

    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for ( float t = minT + step; t < maxT; t += step ) {
        auto cosOfT = cosf(t);
        float x = 2.0f * ( a * cosOfT + b ) * cosOfT;
        float y = 2.0f * ( a * cosOfT + b ) * sinf (t);

        glVertex2f(x - 0.3f, y);
    }

    glEnd();
}

void Display5() {
    constexpr float a = 0.2f;

    constexpr float step = 0.05f;
    constexpr auto pi = static_cast <float> (M_PI);
    constexpr auto minT = - (pi / 2.0f);
    constexpr auto maxT = pi / 2.0f;
    constexpr float excludedT[] = { - (pi / 6.0f), pi / 6.0f };

    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for ( float t = minT + step; t < maxT; t += step ) {
        for ( auto outT : excludedT )
            if ( t == outT )
                continue;

        float sqCosOfT = powf(cosf(t), 2.0f);
        float n = 4 * sqCosOfT - 3;

        float x = a / n;
        float y = (a * tanf(t)) / n;

        glVertex2f(x, y);
    }

    glEnd();
}

void Display6() {
    constexpr float a = 0.1f, b = 0.2f;

    constexpr float step = 0.05f;
    constexpr auto minT = -15.0f;
    constexpr auto maxT = 15.0f;

    glColor3f(1,0.1,0.1); // rosu

    glBegin(GL_LINE_STRIP);

    for ( float t = minT; t < maxT; t += step ) {
        float x = a * t - b * sinf(t);
        float y = a - b * cosf(t);

        glVertex2f(x, y);
    }

    glEnd();
}

void Display7() {
    constexpr float r = 0.3f, R = 0.1f;

    constexpr float step = 0.05f;
    constexpr float minT = 0.0f;
    constexpr float maxT = static_cast<float>(M_PI) * 2.0f;

    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for ( float t = minT; t <= maxT; t += step ) {
        float x = (R + r) * cosf ( r / R * t ) - r * cos ( t + r / R * t );
        float y = (R + r) * sinf ( r / R * t ) - r * sin ( t + r / R * t );

        glVertex2f(x, y);
    }

    glEnd();
}

void Display8() {
    constexpr float r = 0.3f, R = 0.1f;

    constexpr float step = 0.05f;
    constexpr float minT = 0.0f;
    constexpr float maxT = static_cast<float>(M_PI) * 2.0f;

    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for ( float t = minT; t <= maxT; t += step ) {
        float x = (R - r) * cosf ( r / R * t ) - r * cos ( t - r / R * t );
        float y = (R - r) * sinf ( r / R * t ) - r * sin ( t - r / R * t );

        glVertex2f(x, y);
    }

    glEnd();
}

void Display9() {
    constexpr float a = 0.4f;
    constexpr float step = 0.005f;

    constexpr auto pi = static_cast < float > (M_PI);
    constexpr float minT = - ( pi / 4.0f );
    constexpr float maxT = ( pi / 4.0f );

    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for ( float t = minT + step; t < maxT; t += step ) {
        float r = a * sqrtf(2 * cosf(2 * t));

        float x = r * cosf(t);
        float y = r * sinf(t);
        glVertex2f(x, y);
    }

//    glEnd();
//    glBegin(GL_LINE_STRIP);
    for ( float t = maxT - step; t > minT; t -= step ) {
        float r = - a * sqrtf(2 * cosf(2 * t));

        float x = r * cosf(t);
        float y = r * sinf(t);
        glVertex2f(x, y);
    }

    glEnd();
}

void Display0() {
    constexpr float a = 0.02f;
    constexpr float step = 0.05f;

    constexpr float maxT = 100.0f;

    constexpr auto e = static_cast<float>(M_E);

    glColor3f(1,0.1,0.1); // rosu
    glBegin(GL_LINE_STRIP);

    for ( float t = step; t < maxT; t+= step ) {
        float r = a * powf(e, 1.0f + t);
        float x = r * cosf(t);
        float y = r * sinf(t);

        glVertex2f(x, y);
    }

    glEnd();
}

void Init(void) {

    glClearColor(1.0,1.0,1.0,1.0);

    glLineWidth(1);

//   glPointSize(4);

    glPolygonMode(GL_FRONT, GL_LINE);
}

void Display(void) {
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
        case '9':
            Display9();
            break;
        case '0':
            Display0();
            break;
        default:
            break;
    }

    glFlush();
}

void Reshape(int w, int h) {
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}

void KeyboardFunc(unsigned char key, int x, int y) {
    prevKey = key;
    if (key == 27) // escape
        exit(0);
    glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) {
}

int main(int argc, char** argv) {

    glutInit(&argc, argv);

    glutInitWindowSize(dim, dim);

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