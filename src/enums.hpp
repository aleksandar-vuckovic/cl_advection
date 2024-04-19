#ifndef ENUMS_HPP
#define ENUMS_HPP

enum class InitShape
{
    sphere,
    plane,
    paraboloid,
    ellipsoid
};

enum class InitSphereMethod
{
    signedDistance,
    legacySignedSquaredDistance,
    legacySignedScaledSquaredDistance,
};

#endif // ENUMS_HPP
