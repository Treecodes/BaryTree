#include <stdlib.h>
#include <math.h>
#include <stdio.h>

int main()
{

        /* runtime parameters */
        double ox, oy, oz;
        double dx, dy, dz;
        int nx, ny, nz;

        /* input and output files */
        char sampout[20];
        FILE *fp;

        //local variables (these decs may be unnecessary in C)
        int i, j, k;
        double xx, yy, zz;


        printf("Enter name of output file: \n");
        scanf("%s", sampout);
        //sampout = "out.txt";

        printf("Enter origin x, y, z: \n");
        scanf("%lf %lf %lf", &ox, &oy, &oz);
        //ox = -24;
        //oy = ox;
        //oz = ox;

        printf("Enter grid spacing x, y, z: \n");
        scanf("%lf %lf %lf", &dx, &dy, &dz);
        //dx = 0.5;
        //dy = dx;
        //dz = dx;

        printf("Enter number gridpoints nx, ny, nz: \n");
        scanf("%d %d %d", &nx, &ny, &nz);
        //nx = 100;
        //ny = nx;
        //nz = nx/


        /* Writing coordinates for the targets */
        fp = fopen(sampout, "wb");
        if (!fp) {
                perror("File opening failed");
                return EXIT_FAILURE;
        }

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nz; k++)
                {
                    xx = ox + i*dx;
                    yy = oy + j*dy;
                    zz = oz + k*dz;
                    fwrite(&xx, sizeof(double), 1, fp);
                    fwrite(&yy, sizeof(double), 1, fp);
                    fwrite(&zz, sizeof(double), 1, fp);
                }
            }
        }

        fclose(fp);

}
