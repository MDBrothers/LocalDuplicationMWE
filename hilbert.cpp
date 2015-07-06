//http://stackoverflow.com/questions/14519267/algorithm-for-generating-a-3d-hilbert-space-filling-curve-in-python
//answered Nov 22 '14 at 16:59
//by
//kylefinn

function hilbertC(s, x, y, z, dx, dy, dz, dx2, dy2, dz2, dx3, dy3, dz3)
    {
        if(s==1)
        {
            red[m] = x;
            green[m] = y;
            blue[m] = z;
            m++;
        }
        else
        {
            s/=2;
            if(dx<0) x-=s*dx;
            if(dy<0) y-=s*dy;
            if(dz<0) z-=s*dz;
            if(dx2<0) x-=s*dx2;
            if(dy2<0) y-=s*dy2;
            if(dz2<0) z-=s*dz2;
            if(dx3<0) x-=s*dx3;
            if(dy3<0) y-=s*dy3;
            if(dz3<0) z-=s*dz3;
            hilbertC(s, x, y, z, dx2, dy2, dz2, dx3, dy3, dz3, dx, dy, dz);
            hilbertC(s, x+s*dx, y+s*dy, z+s*dz, dx3, dy3, dz3, dx, dy, dz, dx2, dy2, dz2);
            hilbertC(s, x+s*dx+s*dx2, y+s*dy+s*dy2, z+s*dz+s*dz2, dx3, dy3, dz3, dx, dy, dz, dx2, dy2, dz2);
            hilbertC(s, x+s*dx2, y+s*dy2, z+s*dz2, -dx, -dy, -dz, -dx2, -dy2, -dz2, dx3, dy3, dz3);
            hilbertC(s, x+s*dx2+s*dx3, y+s*dy2+s*dy3, z+s*dz2+s*dz3, -dx, -dy, -dz, -dx2, -dy2, -dz2, dx3, dy3, dz3);
            hilbertC(s, x+s*dx+s*dx2+s*dx3, y+s*dy+s*dy2+s*dy3, z+s*dz+s*dz2+s*dz3, -dx3, -dy3, -dz3, dx, dy, dz, -dx2, -dy2, -dz2);
            hilbertC(s, x+s*dx+s*dx3, y+s*dy+s*dy3, z+s*dz+s*dz3, -dx3, -dy3, -dz3, dx, dy, dz, -dx2, -dy2, -dz2);
            hilbertC(s, x+s*dx3, y+s*dy3, z+s*dz3, dx2, dy2, dz2, -dx3, -dy3, -dz3, -dx, -dy, -dz);
        }
    }
    m=0;
    hilbertC(256,0,0,0,1,0,0,0,1,0,0,0,1);
