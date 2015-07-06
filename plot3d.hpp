#ifndef PLOT3D_H
#define PLOT3D_H

#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include<X11/Xlib.h>
#include<X11/XKBlib.h>
#include<GL/glx.h>
#include<GL/glext.h>
#include<GL/glu.h>
#include<iostream>
#include<sstream>
#include<string>
#include<exception>
#include<algorithm>
#include <cstdlib>

int plotWithNeighborhoods(const std::vector<float>& vertices, std::vector<std::vector<int> >& neighborhoods);
#endif
