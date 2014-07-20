#ifndef MATRIXHELPER_HH
#define MATRIXHELPER_HH

std::list<std::pair< AlignmentFunctor , AlignmentFunctor > > couplings;


template< class Couplings, class LocalMatrix>
void axpyCouplings( const Couplings& couplings,
                    const LocalMatrix& localMatrix,
                    const size_t local_i,
                    const size_t local_j)
{
  for( auto coupling : couplings)
    { 
      auto i = coupling.first(local_i);
      auto j = coupling.second(local_j);

      localMatrix( i,j).add( fluxValue[ i ][ j ]*phi[ i ]*phi[ j ]*factor);
    }
}
class AlignmentFunctor
{

public:
  AlignmentFunctor( size_t component, size_t dimension):
  component_(component),
  dimension_(dimension){}

  int operator()( size_t local_index)
  {
    return local_index*dimension_*component_;
  }

private:
  int component_;
  int dimensont_;
};






#endif
