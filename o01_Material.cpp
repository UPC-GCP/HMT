// Imports
#include <json/json.h>

// Self-Imports
#include "o01_Material.h"

Material::Material(Json::Value materials, double gravity){
    // List
    vMat.resize(materials.size());

    // Store Materials
    for (Json::Value::ArrayIndex i = 0; i < materials.size(); i++){
        vMat[i].rho = materials[i]["rho"].asDouble();
		vMat[i].gamma = materials[i]["gamma"].asDouble();
        vMat[i].cp = materials[i]["cp"].asDouble();
        vMat[i].mu = materials[i]["mu"].asDouble();
        vMat[i].beta = materials[i]["beta"].asDouble();
    }

    // External Properties
    g = gravity;
}

void Material::setInitialConditions(double initPhi, Json::Value initVF){ // PHISolver (Value)
    // Main Initial Conditions
    Phi0 = initPhi;
    
    // Velocity Field Initial Conditions
    sVF0.resize(initVF.size());
    for (Json::Value::ArrayIndex i = 0; i < initVF.size(); i++) {sVF0[i] = initVF[i].asString();}
}

void Material::setInitialConditions(std::string pathPhi, Json::Value initVF){ // PHISolver (Path)
    // Solver Variable
    sPhi0 = pathPhi; bPath = true;

    // Velocity Field
    sVF0.resize(initVF.size());
    for (Json::Value::ArrayIndex i = 0; i < initVF.size(); i++) {sVF0[i] = initVF[i].asString();}
}

void Material::setInitialConditions(double initT, double initP, Json::Value initVF){ // NSSolver (Value)
    // Solver Variable
    T0 = initT; P0 = initP;

    // Velocity Field
    VF0.resize(initVF.size());
    for (Json::Value::ArrayIndex i = 0; i < initVF.size(); i++) {VF0[i] = initVF[i].asDouble();}
}

void Material::setInitialConditions(std::string pathT, std::string pathP, Json::Value pathVF){ // NSSolver (Path)
    // Solver Variable
    sT0 = pathT; sP0 = pathP; bPath = true;

    // Velocity Field
    sVF0.resize(pathVF.size());
    for (Json::Value::ArrayIndex i = 0; i < pathVF.size(); i++) {sVF0[i] = pathVF[i].asString();}
}

void Material::setInitialConditions(double initPhi){ // Old: PHISolver 
    // Solver Variable
    Phi0 = initPhi;
}
