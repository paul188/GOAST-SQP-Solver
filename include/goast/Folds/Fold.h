#include <vector>
#include <tinysplinecxx.h>

class Fold {
    public: 
        Fold();
        ~Fold();
    private:
        // spline to be evaluated and adapted to the mesh
        tsBSpline spline;
};

class ExplicitFold : public Fold {
    public:
        ExplicitFold();
        ~ExplicitFold();
};

class ImplicitFold : public Fold {
    public:
        ImplicitFold();
        ~ImplicitFold();
};