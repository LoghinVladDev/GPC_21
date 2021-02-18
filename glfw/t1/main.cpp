//
// Created by loghin on 18.02.2021.
//
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <cstdio>
#include <cstdlib>

//#include <Shader.h>

#define _W_WIDTH 800
#define _W_HEIGHT 600

auto keyboardInputCallback (GLFWwindow *, int, int, int, int) noexcept -> void;
auto mouseMoveCallback (GLFWwindow *, double, double) noexcept -> void;
auto mouseInputCallback (GLFWwindow *, int, int, int) noexcept -> void;
auto mouseScrollCallback (GLFWwindow *, double, double) noexcept -> void;
auto mouseEnterLeaveCallback (GLFWwindow *, int) noexcept -> void;
auto framebufferResizeCallback (GLFWwindow *, int, int) noexcept -> void;
auto errorCallback (int, const char *) -> void;

auto update(double) -> void;
auto render(double) -> void;

constexpr const char * vertexShaderSrc =
        "#version 330 core\n"
        "layout (location = 0) in vec3 vertexPos;\n"
        "void main() {\n"
        "gl_Position = vec4(vertexPos.x, vertexPos.y, vertexPos.z, 1.0);\n"
        "}\0";

constexpr const char * fragmentShaderSrc =
        "#version 330 core\n"
        "out vec4 fragColor;"
        "void main() {\n"
        "fragColor = vec4(";

unsigned int vertexShader;
unsigned int fragmentShader;

void initShaders();

int lastPressedKey;

int main (int argc, char ** argv) {
    if ( GLFW_TRUE != glfwInit() ) {
        printf("Failed to initialise GLFW\n");
        return 0;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow * window = glfwCreateWindow(_W_WIDTH, _W_HEIGHT, argv[0], nullptr, nullptr);

    if ( window == nullptr ) {
        printf("Failed to create GLFW window\n");
        glfwTerminate();
        return 0;
    }

    glfwMakeContextCurrent(window);

    if ( ! gladLoadGLLoader((GLADloadproc) glfwGetProcAddress )) {
        printf("Failed to initialise GLAD\n");
        glfwTerminate();
        return 0;
    }

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    glfwSetCursorPosCallback(window, mouseMoveCallback);
    glfwSetMouseButtonCallback(window, mouseInputCallback);
    glfwSetScrollCallback(window, mouseScrollCallback);
    glfwSetCursorEnterCallback(window, mouseEnterLeaveCallback);
    glfwSetErrorCallback(errorCallback);
    glfwSetFramebufferSizeCallback(window, framebufferResizeCallback);
    glfwSetKeyCallback(window, keyboardInputCallback);

    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glLineWidth(3);
    glPointSize(4);

    glPolygonMode(GL_FRONT, GL_LINE);

    initShaders();

    double lastFrame = 0.0f;
    while ( ! glfwWindowShouldClose(window)) {
        double currentFrame = glfwGetTime();
        double deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        glClear(GL_COLOR_BUFFER_BIT);

        update(deltaTime);
        render(deltaTime);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

auto keyboardInputCallback (
        GLFWwindow * pWindow,
        int key,
        int scanCode,
        int action,
        int modifierFlags
) noexcept -> void {
    if ( action == GLFW_REPEAT )
        return;

    if ( action == GLFW_PRESS ) {
        lastPressedKey = key;


        if ( key == GLFW_KEY_ESCAPE )
            glfwSetWindowShouldClose(pWindow, true);

    } else if ( action == GLFW_RELEASE ) {

    }

    printf("[KB] You %s : <%c>. AdvData : k = %d; c = %d; a = %d; mF = %d\n",
        (action == GLFW_PRESS ? "pressed" : "released" ),
        static_cast<unsigned char>(key),
        key,
        scanCode,
        action,
        modifierFlags
    );
}

auto mouseMoveCallback (GLFWwindow * pWindow, double x, double y) noexcept -> void {
//    printf("[M] You moved mouse to (%lf, %lf)\n",
//        x,
//        y
//    );
}

auto mouseInputCallback (GLFWwindow * pWindow, int button, int action, int modifierFlags) noexcept -> void {
    printf("[M] You %s button %s. AdvData : b = %d; a = %d; mF = %d\n",
           ( action == GLFW_PRESS ? "pressed" : "released" ),
           ( button == GLFW_MOUSE_BUTTON_LEFT ? "left" :
                ( button == GLFW_MOUSE_BUTTON_RIGHT ? "right" :
                  ( button == GLFW_MOUSE_BUTTON_MIDDLE ? "middle" : "other" )
                )
           ),
           button, action, modifierFlags
    );
}

auto mouseScrollCallback ( GLFWwindow * pWindow, double xOffset, double yOffset ) noexcept -> void {
    printf("[M] You scrolled %lf\n", yOffset);
}

auto mouseEnterLeaveCallback ( GLFWwindow * pWindow, int entered ) noexcept -> void {
    if ( entered ) {
        printf("Mouse Entered Window\n");
    } else {
        printf("Mouse Left Window\n");
    }
}

auto framebufferResizeCallback (GLFWwindow * pWindow, int newW, int newH) noexcept -> void {
    printf("[WReshape] Resized Window to (%d, %d)\n", newW, newH);
    glViewport(0, 0, newW, newH);
}


auto errorCallback (int err, const char * text) -> void {
    printf("[ECallback] Error Encountered : %d. %s\n", err, text);
}

void Display1();
void Display2();
void Display3();
void Display4();
void Display5();
void Display6();
void Display7();
void Display8();

void (* displayLogicFunction) ();

auto update(double) -> void {
#define SWITCH_DISPLAY_LOGIC(_key) case _key: displayLogicFunction = Display ## _key; break;
    switch(lastPressedKey - '0') {
        SWITCH_DISPLAY_LOGIC(1)
        SWITCH_DISPLAY_LOGIC(2)
        SWITCH_DISPLAY_LOGIC(3)
        SWITCH_DISPLAY_LOGIC(4)
        SWITCH_DISPLAY_LOGIC(5)
        SWITCH_DISPLAY_LOGIC(6)
        SWITCH_DISPLAY_LOGIC(7)
        SWITCH_DISPLAY_LOGIC(8)
        default: displayLogicFunction = nullptr;
    }
}

auto render(double) -> void {
    if ( displayLogicFunction != nullptr )
        displayLogicFunction ();
}


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
    glColor3f(1,0.1,0.1);
    glBegin(GL_POINTS);
    glVertex2f(0.5f, 0.5f);
    glVertex2f(0.5f, -0.5f);
    glVertex2f(-0.5f, 0.5f);
    glVertex2f(-0.5f, -0.5f);
    glEnd();
}

void Display4() {
    glColor3f(1,0.1,0.1);
    glBegin(GL_LINE_STRIP);
    glVertex2f(1.0f, 1.0f);
    glVertex2f(1.0f, 0.85f);
    glVertex2f(0.6f, 0.65f);
    glVertex2f(0.6f, 0.5f);
    glEnd();
}

void Display5() {
    glColor3f(1,0.1,0.1);
    glBegin(GL_LINE_LOOP);
    glVertex2f(1.0f, 1.0f);
    glVertex2f(1.0f, 0.85f);
    glVertex2f(0.6f, 0.65f);
    glVertex2f(0.6f, 0.5f);
    glEnd();
}

void Display6() {
    glPolygonMode(GL_FRONT, GL_FILL);
    glColor3f(1,0.1,0.1);
    glBegin(GL_TRIANGLES);
    glVertex2f(-1.0f, -1.0f);
    glVertex2f(-0.8f, -0.8f);
    glVertex2f(-1.0f, -0.8f);

    glVertex2f(0.8f, 0.8f);
    glVertex2f(1.0f, 0.8f);
    glVertex2f(1.0f, 1.0f);
    glEnd();
}

void Display7() {
    glColor3f(1,0.1,0.1);
    glPolygonMode(GL_FRONT, GL_FILL);
    glBegin(GL_QUADS);
    glVertex2f(1.0f, 1.0f);
    glVertex2f(0.3f, 0.7f);
    glVertex2f(0.5f, 0.5f);
    glVertex2f(1.0f, 0.5f);
    glEnd();
}

#include <cmath>
void Display8() {
    struct Point { float x; float y; Point(float x, float y) : x(x), y(y) {} };
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

    drawHex(0.2f, 0.15f, 0.88f, 0.75f);
    drawHex(1.0f, 0.1f, 0.1f, 0.5f);
    drawHex(1.0f, 1.0f, 1.0f, 0.48f);
}

void initShaders() {
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSrc, nullptr);
    glCompileShader(vertexShader);

    int success;
    char error[512];

    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);

    if ( ! success ) {
        glGetShaderInfoLog ( vertexShader, 512, nullptr, error );
        printf("Vertex Shader Compilation Error : %s\n", error);
    }
}