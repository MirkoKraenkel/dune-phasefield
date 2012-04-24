#ifndef LIMITER_DGTESTLOOP_HH
#define LIMITER_DGTESTLOOP_HH

#ifdef ONLY_ONE_P
#if POLORDER == 0 || POLORDER == -1
#define DG_ONE_P DG_P0
#elif POLORDER == 1 
#define DG_ONE_P DG_P1
#elif POLORDER == 2 
#define DG_ONE_P DG_P2
#elif POLORDER == 3 
#define DG_ONE_P DG_P3
#elif POLORDER == 4 
#define DG_ONE_P DG_P4
#elif POLORDER == 5 
#define DG_ONE_P DG_P5
#elif POLORDER == 6 
#define DG_ONE_P DG_P6
#elif POLORDER == 7 
#define DG_ONE_P DG_P7
#elif POLORDER == 8 
#define DG_ONE_P DG_P8
#endif
namespace DG_ONE_P {
  void simulate();
}
#else // for several polynomial orders 
namespace DG_P0 {
  void simulate();
}
namespace DG_P1 {
  void simulate();
}
namespace DG_P2 {
  void simulate();
}
namespace DG_P3 {
  void simulate();
}
namespace DG_P4 {
  void simulate();
}
#endif

#endif
