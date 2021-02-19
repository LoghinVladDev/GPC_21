#include <GL/glut.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>

// dimensiunea ferestrei in pixeli
#define dim 300
// numarul maxim de iteratii pentru testarea apartenentei la mult.Julia-Fatou
#define NRITER_JF 5000
// modulul maxim pentru testarea apartenentei la mult.Julia-Fatou
#define MODMAX_JF 10000000
// ratii ptr. CJuliaFatou
#define RX_JF 0.01
#define RY_JF 0.01

unsigned char prevKey;

class CComplex {
public:
    CComplex() : re(0.0), im(0.0) {}
    CComplex(double re1, double im1) : re(re1 * 1.0), im(im1 * 1.0) {}
    CComplex(const CComplex &c) : re(c.re), im(c.im) {}
    ~CComplex() {}

    CComplex &operator=(CComplex const&c)
    {
        re = c.re;
        im = c.im;
        return *this;
    }

    double getRe() {return re;}
    void setRe(double re1) {re = re1;}

    double getIm() {return im;}
    void setIm(double im1) {im = im1;}

    double getModul() {return sqrt(re * re + im * im);}

    int operator==(CComplex &c1)
    {
        return ((re == c1.re) && (im == c1.im));
    }

    CComplex pow2()
    {
        CComplex rez;
        rez.re = powl(re * 1.0, 2) - powl(im * 1.0, 2);
        rez.im = 2.0 * re * im;
        return rez;
    }

    friend CComplex operator+(CComplex const&c1, CComplex const&c2);
    friend CComplex operator*(CComplex const&c1, CComplex const&c2);

    void print(FILE *f)
    {
        fprintf(f, "%.20f%+.20f i", re, im);
    }

private:
    double re, im;
};

CComplex operator+(CComplex const &c1, CComplex const &c2)
{
    CComplex rez(c1.re + c2.re, c1.im + c2.im);
    return rez;
}

CComplex operator*(CComplex const&c1, CComplex const&c2)
{
    CComplex rez(c1.re * c2.re - c1.im * c2.im,
                 c1.re * c2.im + c1.im * c2.re);
    return rez;
}

class CJuliaFatou {
public:
    CJuliaFatou()
    {
        // m.c se initializeaza implicit cu 0+0i

        m.nriter = NRITER_JF;
        m.modmax = MODMAX_JF;
    }

    CJuliaFatou(CComplex &c)
    {
        m.c = c;
        m.nriter = NRITER_JF;
        m.modmax = MODMAX_JF;
    }

    ~CJuliaFatou() {}

    void setmodmax(double v) {assert(v <= MODMAX_JF); m.modmax = v;}
    double getmodmax() {return m.modmax;}

    void setnriter(int v) {assert(v <= NRITER_JF); m.nriter = v;}
    int getnriter() {return m.nriter;}

    // testeaza daca x apartine multimii Julia-Fatou Jc
    // returneaza 0 daca apartine, -1 daca converge finit, +1 daca converge infinit
    int isIn(CComplex &x)
    {
        int rez = 0;
        // tablou in care vor fi memorate valorile procesului iterativ z_n+1 = z_n * z_n + c
        CComplex z0,z1;

        z0 = x;
        for (int i = 1; i < m.nriter; i++)
        {
            z1 = z0 * z0 + m.c;
            if (z1 == z0)
            {
                // x nu apartine m.J-F deoarece procesul iterativ converge finit
                rez = -1;
                break;
            }
            else if (z1.getModul() > m.modmax)
            {
                // x nu apartine m.J-F deoarece procesul iterativ converge la infinit
                rez = 1;
                break;
            }
            z0 = z1;
        }

        return rez;
    }

    // afisarea multimii J-F care intersecteaza multimea argument
    void display(double xmin, double ymin, double xmax, double ymax)
    {
        glPushMatrix();
        glLoadIdentity();

//    glTranslated((xmin + xmax) * 1.0 / (xmin - xmax), (ymin + ymax)  * 1.0 / (ymin - ymax), 0);
//    glScaled(1.0 / (xmax - xmin), 1.0 / (ymax - ymin), 1);
        // afisarea propriu-zisa
        glBegin(GL_POINTS);
        for(double x = xmin; x <= xmax; x+=RX_JF)
            for(double y = ymin; y <= ymax; y+=RY_JF)
            {
                CComplex z(x, y);
                int r = isIn(z);
//        z.print(stdout);
                if (r == 0)
                {
//          fprintf(stdout, "   \n");
                    glVertex3d(x,y,0);
                }
                else if (r == -1)
                {
//          fprintf(stdout, "   converge finit\n");
                }
                else if (r == 1)
                {
//          fprintf(stdout, "   converge infinit\n");
                }
            }
        fprintf(stdout, "STOP\n");
        glEnd();

        glPopMatrix();
    }

private:
    struct SDate {
        CComplex c;
        // nr. de iteratii
        int nriter;
        // modulul maxim
        double modmax;
    } m;
};

// multimea Julia-Fatou pentru z0 = 0 si c = -0.12375+0.056805i
void Display1() {
    CComplex c(-0.12375, 0.056805);
    CJuliaFatou cjf(c);

    glColor3f(1.0, 0.1, 0.1);
    cjf.setnriter(30);
    cjf.display(-0.8, -0.4, 0.8, 0.4);
}

// multimea Julia-Fatou pentru z0 = 0 si c = -0.012+0.74i
void Display2() {
    CComplex c(-0.012, 0.74);
    CJuliaFatou cjf(c);

    glColor3f(1.0, 0.1, 0.1);
    cjf.setnriter(30);
    cjf.display(-1, -1, 1, 1);
}

#include <utility>
class Complex {
private:
    float _re {0.0f};
    float _im {0.0f};

public:
    explicit constexpr Complex (float re = 0.0f, float im = 0.0f) noexcept : _re(re), _im(im) {}
    constexpr Complex ( Complex const & o ) noexcept = default;
    constexpr Complex ( Complex && o ) noexcept : _re(std::exchange(o._re, 0.0f)), _im(std::exchange(o._im, 0.0f)) {}

    constexpr Complex & operator = ( Complex const & o ) noexcept {
        if ( this == & o ) return * this;

        this->_re = o._re;
        this->_im = o._im;

        return * this;
    }

    constexpr auto operator + ( Complex const & o ) const noexcept -> Complex {
        return Complex ( this->_re + o._re, this->_im + o._im );
    }

    constexpr auto operator - ( Complex const & o ) const noexcept -> Complex {
        return Complex ( this->_re - o._re, this->_im - o._im );
    }

    constexpr auto operator * ( Complex const & o ) const noexcept -> Complex {
        return Complex (
                this->_re * o._re - this->_im * o._im,
                this->_re * o._im + this->_im * o._re
        );
    }

    constexpr auto operator / ( Complex const & o ) const noexcept -> Complex {
        auto divisor = o._re * o._re + o._im * o._im;
        return Complex (
                (this->_re * o._re + this->_im * o._im) / divisor,
                (this->_im * o._re - this->_re * o._im) / divisor
        );
    }

    constexpr auto mod () const noexcept -> float {
        return sqrtf( this->_re * this->_re + this->_im * this->_im );
    }

    constexpr auto conjugate () const noexcept -> Complex {
        return Complex ( this->_re, - this->_im );
    }

    constexpr auto operator == ( Complex const & o ) const noexcept -> bool {
        return this->_re == o._re && this->_im == o._im;
    }

    constexpr auto operator != ( Complex const & o ) const noexcept -> bool {
        return ! ( this->_re == o._re && this->_im == o._im );
    }

    constexpr static auto hash( Complex const & o, float clusterPrecision, float additive = 100.0f ) noexcept -> std::size_t {
        return static_cast<std::size_t>(o._re * clusterPrecision * additive + o._im * clusterPrecision);
    }

    constexpr auto getRe() const noexcept -> float {
        return _re;
    }

    constexpr auto getIm() const noexcept -> float {
        return _im;
    }
};

#include <unordered_map>
#include <thread>
#include <vector>

namespace std {
    template <>
    struct hash <Complex> {

        std::size_t operator () ( Complex const & k ) const {
            return Complex::hash(k, 100.0f);
        }
    };
}

#include <atomic>
class Mandelbrot {
private:
    int _iterCount;
    float _precision;
    int _threadCount;

//    std::unordered_map < Complex, std::pair < bool, int > > _set;
    std::unordered_multimap < Complex, std::pair < bool, int > > _set;

public:
    Mandelbrot() noexcept = delete;
    Mandelbrot(Mandelbrot const &) noexcept = default;
    explicit Mandelbrot ( int iterCount, float prec = 0.05f, int threadCnt = 1 ) noexcept : _iterCount(iterCount), _precision(prec), _threadCount(threadCnt) {}

    auto isIn (Complex const & p) const noexcept -> std::pair< bool, int > {
        Complex z;

        for ( int i = 0; i < this->_iterCount; i ++ ) {
            if ( z.mod() > 2.0f ) return { false, i };

            z = z * z + p;
        }

        return {true, -1};
    }

    class SolverThread {
    private:
//        std::atomic <bool> done { true };
//        std::atomic <bool> firstRun { true };
//        std::atomic <bool> start { false };
//        std::atomic <bool> shouldRun { true };
        std::atomic <bool> done;
        std::atomic <bool> firstRun;
        std::atomic <bool> start;
        std::atomic <bool> shouldRun;

        Complex c;
        int iterCount;

        std::pair < bool, int > out;

        void inRun () {
            Complex z;

            for ( int i = 0 ;i < this->iterCount; i++ ) {
                if ( z.mod() > 2.0f ) {
                    this->out = {false, i};
                    this->done = true;
                    return;
                }

                z = z * z + this->c;
            }

            this->out = {true, -1};
            this->done = true;
        }

    public:

        SolverThread(SolverThread const & o) noexcept {
            this->done = static_cast<bool>(o.done);
            this->shouldRun = static_cast<bool>(o.shouldRun);
            this->start = static_cast<bool>(o.start);
            this->firstRun = static_cast<bool>(o.firstRun);
            this->iterCount = o.iterCount;
            this->out = o.out;
            this->c = o.c;
        }

        SolverThread() noexcept {
            this->done = true;
            this->firstRun = true;
            this->start = false;
            this->shouldRun = true;
        }

        void run() {
            while ( this->shouldRun ) {
                while ( ! this->start );

                this->done = false;
                this->inRun();
                this->start = false;
            }
        }

        auto startTask(Complex const & n, int it) noexcept -> void {
            this->iterCount = it;
            this->c = n;

//            this->done = false;
            this->firstRun = false;
            this->start = true;
        }

        auto stop () noexcept -> void {
            this->shouldRun = false;
            this->start = true;
        }

        auto getOut() const -> std::pair < bool, int > const & {
            return out;
        }

        auto getC() const -> Complex const & {
            return c;
        }

        auto isDone() const noexcept -> std::atomic<bool> const & { return this->done; }
        auto isFirstRun () const noexcept -> std::atomic<bool> const & { return this->firstRun; }
    };

    // red - orange - yellow - green - blue
    // 1.0, 0.0, 0.0,
    //// 1,0, 0.5, 0.0
    // 1.0, 1.0, 0.0,
    //// 0.5, 1.0, 0.0,
    // 0.0, 1.0, 0.0,
    //// 0.0, 1.0, 0.5,
    // 0.0, 1.0, 1.0
    //// 0.0, 0.5, 1.0,
    // 0.0, 0.0, 1.0

    // 1.0, 0.0, 0.0,
    // 1.0, 1.0, 0.0,
    // 0.0, 1.0, 0.0,
    // 0.0, 1.0, 1.0
    // 0.0, 0.0, 1.0

    auto compute () noexcept -> void {
        if ( this->_threadCount == 1 )
            for ( float i = -2.0f; i < 2.0f; i+= this->_precision )
                for ( float j = -2.0f; j < 2.0f; j+= this->_precision )
                    this->_set.emplace( Complex (i, j), this->isIn( Complex(i, j) ) );
        else if ( this->_threadCount > 1 ) {
            std::vector < SolverThread > solvers;
            solvers.resize(this->_threadCount);
//            for ( int i = 0; i < this->_threadCount; i++ )
//                solvers.emplace_back(SolverThread());

            std::vector < std::thread > threads;
            threads.resize(this->_threadCount);
            for ( int i = 0; auto & t : threads )
                t = std::thread(& SolverThread::run, & solvers[i++]);
//                threads.emplace_back(& SolverThread::run, & solvers[i]);

            for ( float i = -2.0f; i < 2.0f; i+= this->_precision )
                for ( float j = -2.0f; j < 2.0f; j+= this->_precision ) {
                    bool solverFound = false;
                    while ( ! solverFound ) {
                        for ( SolverThread & st : solvers ) {
                            if ( st.isDone() ) {
                                if ( ! st.isFirstRun() )
                                    this->_set.emplace( st.getC(), st.getOut() );
                                st.startTask( Complex (i, j), this->_iterCount );
                                solverFound = true;
                                break;
                            }
                        }
                    }
                }

            for ( SolverThread & st : solvers ) {
                while ( ! st.isDone() );
                this->_set.emplace( st.getC(), st.getOut() );
                st.stop();
            }

            for ( std::thread & t : threads )
                t.join();

        }
    }

    auto requiresCompute() const noexcept -> bool {
        return this->_set.empty();
    }

    int getIterCount() const {
        return _iterCount;
    }

    inline auto getData () const noexcept -> std::unordered_multimap < Complex, std::pair < bool, int > > { return this->_set; }
};

#include <map>

int iterCount = 80;
int prevIterCount = iterCount;

Mandelbrot mb(iterCount, 0.003f, 6);

void Display3() {
    if ( iterCount != prevIterCount || mb.requiresCompute() ) {
        prevIterCount = iterCount;
        mb.compute();
    }

    auto clampToViewport = [] (Complex const & c) -> Complex { return Complex(c.getRe() / 2.0f, c.getIm() / 2.0f); };

    glPointSize(1.0f);
    glBegin(GL_POINTS);

    std::map < int, std::tuple < float, float, float > > colorMapping;

    float dist = 1.0f / ( static_cast<float>(mb.getIterCount()) / 4.0f );

    for ( int i = 0; i < mb.getIterCount() / 4; i++ ) colorMapping.emplace( i, std::tuple<float, float, float>{ 1.0f, static_cast<float>(i) * dist , 0.0f } );
    for ( int i =  mb.getIterCount() / 4; i < mb.getIterCount() / 2; i++ ) colorMapping.emplace( i, std::tuple<float, float, float>{ 1.0f - static_cast<float>(i - mb.getIterCount() / 4) * dist, 1.0f, 0.0f } );
    for ( int i = mb.getIterCount() / 2; i < mb.getIterCount() / 4 * 3; i++ ) colorMapping.emplace( i, std::tuple<float, float, float>{ 0.0f, 1.0f, static_cast<float>(i - mb.getIterCount() / 2) * dist } );
    for ( int i = mb.getIterCount() / 4 * 3; i < mb.getIterCount(); i++ ) colorMapping.emplace ( i, std::tuple<float, float, float>{ 0.0f, 1.0f - static_cast<float>(i - mb.getIterCount() / 4 * 3) * dist, 1.0f } );

    for ( auto & e : mb.getData() ) {

        if ( e.second.first ) {
            glColor3f ( 0.0f, 0.0f, 0.0f );
        } else {
            glColor3f (
                std::get<0>(colorMapping.find(e.second.second)->second),
                std::get<1>(colorMapping.find(e.second.second)->second),
                std::get<2>(colorMapping.find(e.second.second)->second)
            );
        }

        auto clamped = clampToViewport(e.first);
        glVertex2f( clamped.getRe(), clamped.getIm() );
    }

    glEnd();
}

struct Point {
    float x, y;

    constexpr Point ( float x, float y ) noexcept: x(x), y(y) {}
    Point ( Point const & p ) noexcept = default;

    constexpr auto operator == (Point const & b) const noexcept -> bool {
        return this->x == b.x && this->y == b.y;
    }

    auto static distance ( Point const & a, Point const & b ) noexcept -> float {
        float xDif = a.x - b.x;
        float yDif = a.y - b.y;
        return sqrtf( xDif * xDif + yDif * yDif );
    }

    auto static angleToHAxisRad ( Point const & a, Point const & b ) noexcept -> float {
        return atan2f(b.y - a.y, b.x - a.x);
    }

    auto inline static angleToHAxisDeg ( Point const & a, Point const & b ) noexcept -> float {
        constexpr static float rtp = 180.0f / static_cast<float> (M_PI);
        return Point::angleToHAxisRad(a, b) * rtp;
    }

    auto inline distance ( Point const & a ) const noexcept -> float { return Point::distance( * this, a ); }
    auto inline angleToHAxisRad( Point const & a ) const noexcept -> float { return Point::angleToHAxisRad( *this, a ); }
    auto inline angleToHAxisDeg( Point const & a ) const noexcept -> float { return Point::angleToHAxisDeg( *this, a ); }
};

namespace std {
    template <>
    struct hash <Point> {
        std::size_t operator () ( Point const & k ) const {
            return static_cast<std::size_t>(k.x * 1000.0f + k.y * 100.0f);
        }
    };
}

class DrawableFractal {
protected:
    std::vector < Point > _points;
    int prevIter;
    int actualIter;
public:
    virtual auto requiresCompute () noexcept -> bool {
        return this->_points.empty() || prevIter != actualIter;
    }

    virtual auto compute (int) noexcept -> void = 0;

    virtual auto draw () noexcept -> void = 0;
};

#include <unordered_set>
class TurtleDrawable {
protected:
    virtual auto nextIteration ( std::vector < Point > const & ) noexcept -> std::vector < Point > = 0;
};

#include <array>
class KochCurve : public DrawableFractal, public TurtleDrawable {
private:
    std::array < Point, 3 > _startPoints;

public:
    KochCurve() noexcept = delete;
    explicit KochCurve( Point const & a, Point const & b, Point const & c ) noexcept : _startPoints ( {a, b, c} ) {}
    explicit KochCurve( std::array < Point, 3 > const & sPoints ) noexcept : _startPoints(sPoints) {}

    auto draw () noexcept -> void final {
        glColor3f(1.0f, 0.1f, 0.1f);
        glBegin(GL_LINE_LOOP);

        for ( auto & p : this->_points )
            glVertex2f( p.x, p.y );

        glEnd();
    }

    auto compute ( int iterCount ) noexcept -> void final {
        this->actualIter = iterCount;
        if ( ! this->requiresCompute() ) {
            return;
        }
        this->prevIter = iterCount;

        this->_points = std::vector < Point > (this->_startPoints.begin(), this->_startPoints.end());
        for ( int i = 0; i < iterCount; i++ ) {
            std::vector < Point > nextPoints;

            for ( int i = 1; i < this->_points.size(); i ++ ) {
                auto subset = this->nextIteration({this->_points[i-1], this->_points[i]});
                for ( auto & p : subset )
                    nextPoints.emplace_back( p );
            }

            auto subset = this->nextIteration({this->_points[this->_points.size() - 1], this->_points[0]});
            for ( auto & p : subset )
                nextPoints.emplace_back( p );

            this->_points = nextPoints;
        }
    }

    auto nextIteration ( std::vector < Point > const & points ) noexcept -> std::vector < Point > final {
        if ( points.size() != 2 )
            return std::vector < Point > ();

        std::vector < Point > nextIterationPoints;

        auto & a = * points.begin();
        auto & b = * (++points.begin());

        float len = a.distance(b) / 3.0f;
        float angle = a.angleToHAxisDeg(b);

        constexpr static auto pi = static_cast<float>(M_PI);
        constexpr auto degToRad = [] (float deg) noexcept -> float { return deg * (pi / 180.0f); };

        Point p1 = { a.x + cosf(degToRad(angle)) * len, a.y + sinf(degToRad(angle)) * len };
        Point p2 = { p1.x + cosf(degToRad(angle - 60.0f)) * len, p1.y + sinf (degToRad(angle - 60.0f)) * len };
        Point p3 = { p2.x + cosf(degToRad(angle + 60.0f)) * len, p2.y + sinf (degToRad(angle + 60.0f)) * len };

        nextIterationPoints.emplace_back(a);
        nextIterationPoints.emplace_back(p1);
        nextIterationPoints.emplace_back(p2);
        nextIterationPoints.emplace_back(p3);
        nextIterationPoints.emplace_back(b);

        return nextIterationPoints;
    }
};

int d4ICount = 1;

void Display4() {
    static KochCurve curve({0.8f, -0.5f}, {0.0f, 0.886f}, {-0.8f, -0.5f});
    curve.compute(d4ICount);
    curve.draw();
}

class DrawableBinaryTreeFractal {
protected:
    struct Node {
        Point p;
        Node * left = nullptr, * right = nullptr;
        float angleWithParent;
        int height;

        ~Node () {
            delete left;
            delete right;
        }
    };

    Node * root = nullptr;
    int prevIter = 0;
    int actualIter = 0;
public:

    virtual auto requiresCompute () noexcept -> bool {
        return this->root == nullptr || prevIter != actualIter;
    }

    virtual auto compute (int) noexcept -> void = 0;
    virtual auto draw () noexcept -> void = 0;

    virtual ~DrawableBinaryTreeFractal() noexcept {
        delete root;
    }
};

class BinaryTree : public DrawableFractal, public TurtleDrawable {

}

void Display5() { }
void Display6() { }
void Display7() { }
void Display8() { }
void Display9() { }
void Display0() { }

void Init(void) {

    glClearColor(1.0,1.0,1.0,1.0);

    glLineWidth(1);

//   glPointSize(3);

    glPolygonMode(GL_FRONT, GL_LINE);
}

char prevDisp = 0;
int d4ICountPrev = 0;

void Display(void) {
    if ( iterCount != prevIterCount && prevDisp == '3' ) {
        glClear(GL_COLOR_BUFFER_BIT);
        Display3();
    }
    if ( prevDisp == '4' ) {
        glClear(GL_COLOR_BUFFER_BIT);
        Display4();
    }

    switch(prevKey) {
        case '1':
            glClear(GL_COLOR_BUFFER_BIT);
            Display1();
            break;
        case '2':
            glClear(GL_COLOR_BUFFER_BIT);
            Display2();
            break;
        case '3':
            glClear(GL_COLOR_BUFFER_BIT);
            prevDisp = '3';
            Display3();
            break;
        case '4':
            glClear(GL_COLOR_BUFFER_BIT);
            Display4();
            prevDisp = '4';
            break;
        case '5':
            glClear(GL_COLOR_BUFFER_BIT);
            Display5();
            break;
        case '6':
            glClear(GL_COLOR_BUFFER_BIT);
            Display6();
            break;
        case '7':
            glClear(GL_COLOR_BUFFER_BIT);
            Display7();
            break;
        case '8':
            glClear(GL_COLOR_BUFFER_BIT);
            Display8();
            break;
        case '9':
            glClear(GL_COLOR_BUFFER_BIT);
            Display9();
            break;
        case '0':
            glClear(GL_COLOR_BUFFER_BIT);
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

    if ( key == 'w' ) {
        iterCount += 5;
        d4ICount++;
    }
    else if ( key == 's' ) {
        if (  iterCount > 5 )
        iterCount -= 5;
        if ( d4ICount > 0 )
            d4ICount --;
    }
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

