#include "plot3d.hpp"

Display                 *dpy;
Window                  root, win;
GLint                   att[]   = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
XVisualInfo             *vi;
GLXContext              glc;
Colormap                cmap;
XSetWindowAttributes    swa;
XWindowAttributes       wa;
XEvent                  xev;
bool terminate = false;

float                   TimeCounter, LastFrameTimeCounter, DT, prevTime = 0.0, FPS;
struct timeval          tv, tv0;
int                     Frame = 1, FramesPerFPS;

GLfloat                 rotation_matrix[16];
float                   rot_z_vel = 0.0, rot_y_vel = 0.0;

class Model {
public:
    int mydimension;
    float bound;
    int arraySize;
    int lengthNeighborhoods;
    const float* mypoints;
    std::vector<std::vector<int> >::iterator neighborhoodsBegin;
    std::vector<std::vector<int> >::iterator neighborhoodsEnd;
    std::vector<float> colorVector;

    float random(float max);

    Model(const std::vector<float>& points, std::vector<std::vector<int> >& neighborhoods);
    void draw();
};

float Model::random(float max) {
    return float(rand())/float(RAND_MAX)*max;
}

Model::Model(const std::vector<float>& points, std::vector<std::vector<int> >& neighborhoods)
{
    mypoints = &points[0];
    arraySize = points.size();
    bound = *std::max_element(points.begin(), points.end()) - *std::min_element(points.begin(), points.end());

    // set color vector
    float c1(Model::random(1.f)), c2(Model::random(1.f)), c3(Model::random(1.f));
    // This is one of many ways to add a color value for root node.
    for(std::vector<std::vector<int> >::iterator hood = neighborhoods.begin(); hood != neighborhoods.end(); ++ hood)
    {
        c1 = Model::random(1.f);
        c2 = Model::random(1.f);
        c3 = Model::random(1.f);
        colorVector.push_back(c1);
        colorVector.push_back(c2);
        colorVector.push_back(c3);
    }

    neighborhoodsBegin = neighborhoods.begin();
    neighborhoodsEnd = neighborhoods.end();
}


void Model::draw() {
    glLineWidth(2.0);
    glPointSize(1.0);
    std::vector<std::vector<int> >::iterator hood = neighborhoodsBegin;

    static int hoodToHighlight = 0;
    glBegin(GL_POINTS);
    glColor3f(1.f,1.f,1.f);
    for(; hood != neighborhoodsEnd; ++ hood) {
        int hoodDistance = std::distance(neighborhoodsBegin, hood);
        for(std::vector<int>::iterator neighbor = hood->begin(); neighbor != hood->end(); ++neighbor) {
            glVertex3f(*(mypoints + *neighbor*3),*( mypoints + *neighbor*3+ 1) ,*( mypoints +  *neighbor*3 + 2));
        }
    }
    glEnd();

    hood = neighborhoodsBegin;

    glBegin(GL_LINES);
    for(; hood != neighborhoodsEnd; ++ hood) {
        if(std::distance(neighborhoodsBegin, hood) == hoodToHighlight) {
            int hoodDistance = std::distance(neighborhoodsBegin, hood);
            glColor3f(colorVector[hoodDistance*3], colorVector[hoodDistance*3 + 1], colorVector[hoodDistance*3 + 2]);

            for(std::vector<int>::iterator neighbor = hood->begin(); neighbor != hood->end(); ++neighbor) {
                glVertex3f(*(mypoints + hoodDistance*3),*( mypoints + hoodDistance*3 + 1) ,*( mypoints +  hoodDistance*3 + 2));
                glVertex3f(*(mypoints + *neighbor*3),*( mypoints + *neighbor*3+ 1) ,*( mypoints +  *neighbor*3 + 2));
            }
        }
    }
    glEnd();

    hoodToHighlight ++;
    hoodToHighlight = hoodToHighlight % std::distance(neighborhoodsBegin, neighborhoodsEnd);
}


void RotateModel() {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glRotatef(rot_y_vel*DT, 0.0, 1.0, 0.0);
    glRotatef(rot_z_vel*DT, 0.0, 0.0, 1.0);
    glMultMatrixf(rotation_matrix);
    glGetFloatv(GL_MODELVIEW_MATRIX, rotation_matrix);
}

void ExposeFunc(Model & myModel) {
    float  aspect_ratio;
    char   info_string[256];

//Resize the viewport
    XGetWindowAttributes(dpy, win, &wa);
    glViewport(0, 0, wa.width, wa.height);
    aspect_ratio = (float)(wa.width) / (float)(wa.height);

//Set up projection (orthographic vs. perspective) and model view (camera location and aim)
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
// glOrtho(-2.50*aspect_ratio, 2.50*aspect_ratio, -2.50, 2.50, -200.0, 200.0);
    gluPerspective(50.0, aspect_ratio, 1.0, 400.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    float eyeX(float(0.0)), eyeY(float(0.0)), eyeZ(-3.0*myModel.bound);
    float aimX(eyeX), aimY(eyeY), aimZ(0.0);
    float upX(0.f), upY(1.f), upZ(0.f);

    gluLookAt(eyeX, eyeY, eyeZ, aimX, aimY, aimZ, upX, upY, upZ);
    glMultMatrixf(rotation_matrix);

//Draw the model
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    myModel.draw();

//Display graphics system stats
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, (float)wa.width, 0, (float)wa.height, -1., 1.);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glColor3f(1.0, 1.0, 1.0);

    sprintf(info_string, "%4.1f seconds * %4.1f fps at %i x %i", TimeCounter, FPS, wa.width, wa.height);
    glRasterPos2i(10, 10);
    glCallLists(strlen(info_string), GL_UNSIGNED_BYTE, info_string);

    sprintf(info_string, "<up,down,left,right> rotate model * <F1> stop rotation ");
    glRasterPos2i(10, wa.height-32);
    glCallLists(strlen(info_string), GL_UNSIGNED_BYTE, info_string);

//Swap display buffer with draw buffer
    glXSwapBuffers(dpy, win);
}

//Create a GL window
void CreateWindow() {

    if((dpy = XOpenDisplay(NULL)) == NULL) {
        printf("\n\tcannot connect to x server\n\n");
        exit(0);
    }

    root = DefaultRootWindow(dpy);

    if((vi = glXChooseVisual(dpy, 0, att)) == NULL) {
        printf("\n\tno matching visual\n\n");
        exit(0);
    }

    if((cmap = XCreateColormap(dpy, root, vi->visual, AllocNone)) == 0) {
        printf("\n\tcannot create colormap\n\n");
        exit(0);
    }

    swa.event_mask = KeyPressMask;
    swa.colormap   = cmap;
    win = XCreateWindow(dpy, root, 0, 0, 700, 700, 0, vi->depth, InputOutput, vi->visual, CWColormap | CWEventMask, &swa);
    XStoreName(dpy, win, "OpenGL Animation");
    XMapWindow(dpy, win);
}

//Setup GL context
void SetupGL() {
    char           font_string[128];
    XFontStruct    *font_struct;

//Create GL context and set as current
    if((glc = glXCreateContext(dpy, vi, NULL, GL_TRUE)) == NULL) {
        printf("\n\tcannot create gl context\n\n");
        exit(0);
    }

    glXMakeCurrent(dpy, win, glc);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.00, 0.00, 0.00, 1.00);

    //Locate a font for text display purposes
    for(int font_size = 14; font_size < 32; font_size += 2) {
        sprintf(font_string, "-adobe-courier-*-r-normal--%i-*", font_size);
        font_struct = XLoadQueryFont(dpy, font_string);

        if(font_struct != NULL) {
            glXUseXFont(font_struct->fid, 32, 192, 32);
            break;
        }
    }

//Initialize the rotation matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glGetFloatv(GL_MODELVIEW_MATRIX, rotation_matrix);
}

//Time counter functions
void InitTimeCounter() {
    gettimeofday(&tv0, NULL);
    FramesPerFPS = 30;
}

void UpdateTimeCounter() {
    LastFrameTimeCounter = TimeCounter;
    gettimeofday(&tv, NULL);
    TimeCounter = (float)(tv.tv_sec-tv0.tv_sec) + 0.000001*((float)(tv.tv_usec-tv0.tv_usec));
    DT = TimeCounter - LastFrameTimeCounter;
}

void CalculateFPS() {
    Frame ++;

    if((Frame%FramesPerFPS) == 0) {
        FPS = ((float)(FramesPerFPS)) / (TimeCounter-prevTime);
        prevTime = TimeCounter;
    }
}

//Cleanup method
void ExitPlotWindow() {
    glXMakeCurrent(dpy, None, NULL);
    glXDestroyContext(dpy, glc);
    XDestroyWindow(dpy, win);
    XCloseDisplay(dpy);
    //exit(0);
}

//User input polling method
void CheckKeyboard() {

    if(XCheckWindowEvent(dpy, win, KeyPressMask, &xev)) {
        char    *key_string = XKeysymToString(XkbKeycodeToKeysym(dpy, xev.xkey.keycode, 0, 0));

        if(strncmp(key_string, "Left", 4) == 0) {
            rot_z_vel -= 0.0*DT;
        }

        else if(strncmp(key_string, "Right", 5) == 0) {
            rot_z_vel += 0.0*DT;
        }
        else if(strncmp(key_string, "Up", 2) == 0) {
            rot_y_vel -= 200.0*DT;
        }

        else if(strncmp(key_string, "Down", 4) == 0) {
            rot_y_vel += 200.0*DT;
        }
        else if(strncmp(key_string, "F1", 2) == 0) {
            rot_y_vel = 0.0;
            rot_z_vel = 0.0;
        }

        else if(strncmp(key_string, "Escape", 5) == 0) {
            ExitPlotWindow();
	    terminate = true;
        }
    }
}

int plotWithNeighborhoods(const std::vector<float>& vertices, std::vector<std::vector<int> >& neighborhoods) {
    Model myModel(vertices, neighborhoods);
    CreateWindow();
    SetupGL();
    InitTimeCounter();

    while(not terminate) {
        UpdateTimeCounter();
        CalculateFPS();
        RotateModel();
        ExposeFunc(myModel);
        usleep(1000);
        CheckKeyboard();
    }

    return 0;
}
