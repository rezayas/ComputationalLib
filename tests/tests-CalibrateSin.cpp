#include "catch.hpp"
#include "iostream"
#include "../include/ComputationalLib/CalibrateSin.hpp"

TEST_CASE("Hello_world_CalibrateSin", "[calibrateSin]") {
    SECTION( "x = 0" ) {
        // x is the initial value
        BoundDouble x{0};

        // x_f will be set to the final outputted value of the model 
        BoundDouble x_f{0};

        F f = [&] (BoundDouble x) {
            return std::sin(x()) + 1; 
        }; 

        x_f = CalibrateSin(x, f);

        double y_f = f(x_f);

        REQUIRE( y_f < 2.001 );
        REQUIRE( y_f > 1.999 );
    }

    SECTION( "x = 1" ) {
        // x is the initial value
        BoundDouble x{1};

        // x_f will be set to the final outputted value of the model 
        BoundDouble x_f{0};

        F f = [&] (BoundDouble x) {
            return std::sin(x()) + 1; 
        }; 

        x_f = CalibrateSin(x, f);

        double y_f = f(x_f);

        REQUIRE( y_f < 2.001 );
        REQUIRE( y_f > 1.999 );
    }
    
    SECTION( "x = 2" ) {
        // x is the initial value
        BoundDouble x{2.0};

        // x_f will be set to the final outputted value of the model 
        BoundDouble x_f{0};

        F f = [&] (BoundDouble x) {
            return std::sin(x()) + 1; 
        }; 

        x_f = CalibrateSin(x, f);

        double y_f = f(x_f);

        REQUIRE( y_f < 2.001 );
        REQUIRE( y_f > 1.999 );
    }
    
    SECTION( "x = 3" ) {
        // x is the initial value
        BoundDouble x{3.0};

        // x_f will be set to the final outputted value of the model 
        BoundDouble x_f{0};

        F f = [&] (BoundDouble x) {
            return std::sin(x()) + 1; 
        }; 

        x_f = CalibrateSin(x, f);

        double y_f = f(x_f);

        REQUIRE( y_f < 2.001 );
        REQUIRE( y_f > 1.999 );
    }
    
    SECTION( "x = 4" ) {
        // x is the initial value
        BoundDouble x{4};

        // x_f will be set to the final outputted value of the model 
        BoundDouble x_f{0};

        F f = [&] (BoundDouble x) {
            return std::sin(x()) + 1; 
        }; 

        x_f = CalibrateSin(x, f);

        double y_f = f(x_f);

        REQUIRE( y_f < 2.001 );
        REQUIRE( y_f > 1.999 );
    }

    SECTION( "x = 5" ) {
        // x is the initial value
        BoundDouble x{5};

        // x_f will be set to the final outputted value of the model 
        BoundDouble x_f{0};

        F f = [&] (BoundDouble x) {
            return std::sin(x()) + 1; 
        }; 

        x_f = CalibrateSin(x, f);

        double y_f = f(x_f);

        REQUIRE( y_f < 1.001 );
        REQUIRE( y_f > 0.999 );
    }

    SECTION( "x = 6" ) {
        // x is the initial value
        BoundDouble x{6};

        // x_f will be set to the final outputted value of the model 
        BoundDouble x_f{0};

        F f = [&] (BoundDouble x) {
            return std::sin(x()) + 1; 
        }; 

        x_f = CalibrateSin(x, f);

        double y_f = f(x_f);

        REQUIRE( y_f < 1.001 );
        REQUIRE( y_f > 0.999 );
    }
}