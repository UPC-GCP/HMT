// Imports
#include <vector> 
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

void Material::setInitialConditions(double initPhi){

    // Main Initial Conditions
    Phi0 = initPhi;
    
}

void Material::setInitialConditions(double initPhi, Json::Value initVF){

    // Main Initial Conditions
    Phi0 = initPhi;
    
    // Velocity Field Initial Conditions
    VF0.resize(initVF.size());
    for (Json::Value::ArrayIndex i = 0; i < initVF.size(); i++){VF0[i] = initVF[i].asDouble();}

}

void Material::setInitialConditions(double initT, double initP, Json::Value initVF){

    // Main Initial Conditions
    T0 = initT; P0 = initP;

    // Velocity Field Initial Conditions
    VF0.resize(initVF.size());
    for (Json::Value::ArrayIndex i = 0; i < initVF.size(); i++){VF0[i] = initVF[i].asDouble();}

}
