#include <GridCreator_NEW.h>

GridCreator_NEW::GridCreator_NEW(InputParser &input_parser,
					    Materials &materials,
					    MPI_Initializer &MPI_communicator):
                        input_parser(input_parser),
                        materials(materials),
                        MPI_communicator(MPI_communicator){
    // Initialize the vectors:
    this->delta_Electromagn[0] = this->input_parser.deltaX;
    this->delta_Electromagn[1] = this->input_parser.deltaY;
    this->delta_Electromagn[2] = this->input_parser.deltaZ;
}