#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <cassert>
#include <cmath>

#define INITIALGRID(I,J,K) this->initialGrid[(I) + this->initialSize[0] * ( (J) + (K) * this->initialSize[1])]

#define NEWGRID(I,J,K)     this->newGrid[    (I) + this->newSize[0]     * ( (J) + (K) * this->newSize[1])]

class ThreeDim_linear_interp{
    private:
        // Initial 3D grid on which we perform trilinear interpolation.
        std::vector<double> &initialGrid;
        std::vector<size_t> &initialSize;
        std::vector<double> &initialStep;
        // New grid:
        std::vector<double> newGrid;
        std::vector<size_t> &newSize;
        std::vector<double> &newStep;
        // Ratio between the two grids:
        size_t ratio;
    public:
        // Contructor:
        ThreeDim_linear_interp(
            std::vector<double> &initialGrid,
            std::vector<size_t> &initialSize,
            std::vector<double> &initialStep,
            std::vector<size_t> &newSize,
            std::vector<double> &newStep,
            size_t ratio) : initialGrid(initialGrid),
                                            initialSize(initialSize),
                                            initialStep(initialStep),
                                            newSize(newSize),
                                            newStep(newStep){
                printf("Welcome to the trilinear grid interpolation...\n");
                this->newGrid.resize(newSize[0]*newSize[1]*newSize[2]);
                this->ratio = ratio;
        }

        // Trilinear interpolation method:
        std::vector<double> trilinearInterp(void){

            printf("Interpolation en cours...\n");

            size_t index;

            size_t use_I = 0;
            size_t use_J = 0;
            size_t use_K = 0;

            double x, y, z;

            double x_0, x_1;
            double y_0, y_1;
            double z_0, z_1;

            double x_d, y_d, z_d;

            double c_000;
            double c_001;
            double c_011;
            double c_010;
            double c_110;
            double c_100;
            double c_111;
            double c_101;

            double c_00;
            double c_01;
            double c_10;
            double c_11;

            double c_0;
            double c_1;
            
            /** LOOP OVER Z DIRECTION **/
            for(size_t K = 0 ; K < this->newSize[2] ; K ++ ){

                if( K != 0 && K%ratio == 0 && K < this->newSize[2] - 1 ){
                    use_K++;
                }

                use_J = 0;

                z_0 = use_K * initialStep[2];
                z_1 = z_0 + initialStep[2];

                z = K * newStep[2];

                z_d = ( z - z_0 ) / ( z_1 - z_0 );

                /** LOOP OVER Y DIRECTION **/
                for(size_t J = 0 ; J < this->newSize[1] ; J ++ ){

                    if( J != 0 && J%ratio == 0 && J < this->newSize[1] - 1 ){
                        use_J++;
                    }

                    use_I = 0;

                    y_0 = use_J * initialStep[1];
                    y_1 = y_0 + initialStep[1];

                    y = J * newStep[1];

                    y_d = ( y - y_0 ) / ( y_1 - y_0 );

                    /** LOOP OVER X DIRECTION **/
                    for(size_t I = 0 ; I < this->newSize[0] ; I ++ ){

                        if( I != 0 && I%ratio == 0 && I < this->newSize[0] - 1 ){
                            use_I++;
                        }

                        x_0 = use_I * initialStep[0];
                        x_1 = x_0 + initialStep[0];

                        x = I * newStep[0];

                        /** COMPUTATION OF TRILINEAR INTERP. COEFFICIENTS **/
                        
                        x_d = ( x - x_0 ) / ( x_1 - x_0 );

                        c_000 = INITIALGRID( use_I   ,use_J   ,use_K   );
                        c_100 = INITIALGRID( use_I+1 ,use_J   ,use_K   );
                        c_110 = INITIALGRID( use_I+1 ,use_J+1 ,use_K   );
                        c_010 = INITIALGRID( use_I   ,use_J+1 ,use_K   );
                        c_011 = INITIALGRID( use_I   ,use_J+1 ,use_K+1 );
                        c_001 = INITIALGRID( use_I   ,use_J   ,use_K+1 );
                        c_111 = INITIALGRID( use_I+1 ,use_J+1 ,use_K+1 );
                        c_101 = INITIALGRID( use_I+1 ,use_J   ,use_K+1 );

                        if( c_000 == 0 ||
                            c_100 == 0 || 
                            c_110 == 0 ||
                            c_010 == 0 ||
                            c_011 == 0 ||
                            c_001 == 0 ||
                            c_111 == 0 ||
                            c_101 == 0){
                                NEWGRID(I,J,K) = 0;
                                continue;
                            }

                        c_00 = c_000 * ( 1 - x_d ) + c_100 * x_d;
                        c_01 = c_001 * ( 1 - x_d ) + c_101 * x_d;
                        c_10 = c_010 * ( 1 - x_d ) + c_110 * x_d;
                        c_11 = c_011 * ( 1 - x_d ) + c_111 * x_d;

                        c_0 = c_00 * ( 1 - y_d ) + c_10 * y_d;
                        c_1 = c_01 * ( 1 - y_d ) + c_11 * y_d;

                        double val = round(c_0 * ( 1 - z_d ) + c_1 * z_d);

                        NEWGRID(I,J,K) = val;

                        if( std::isnan(c_0 * ( 1 - z_d ) + c_1 * z_d) == true ){
                            printf("For (I,J,K) = (%zu,%zu,%zu), I have NaN.\n",I,J,K);
                            std::abort();
                        }

                    }
                }
            }

            return this->newGrid;

        }

};
