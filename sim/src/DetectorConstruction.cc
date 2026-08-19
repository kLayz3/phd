#include "DetectorConstruction.hh"

#include "Config.hh"
#include "G4Box.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

DetectorConstruction::DetectorConstruction()
    : fTargetThickness(cfg::targetArealDensity / cfg::targetDensity),
      fTargetDensity(cfg::targetDensity) {}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    auto* nist = G4NistManager::Instance();
    auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");

    // Isotopically pure 9Be.  For the EM transport in this toy calculation,
    // natural Be would be effectively equivalent, but defining 9Be explicitly
    // makes the intended target unambiguous for later nuclear extensions.
    auto* be9 = new G4Isotope("Be9", 4, 9, 9.0121831 * g / mole);
    auto* beElement = new G4Element("Beryllium9", "Be9", 1);
    beElement->AddIsotope(be9, 100.0 * perCent);
    auto* targetMaterial = new G4Material(
        "Be9Target", fTargetDensity, 1);
    targetMaterial->AddElement(beElement, 1.0);

    const auto worldHalfXY = 25.0 * cm;
    const auto worldHalfZ = 50.0 * cm;
    auto* worldSolid = new G4Box("WorldSolid", worldHalfXY, worldHalfXY, worldHalfZ);
    auto* worldLV = new G4LogicalVolume(worldSolid, vacuum, "WorldLV");
    auto* worldPV = new G4PVPlacement(
        nullptr, {}, worldLV, "WorldPV", nullptr, false, 0, true);

    const auto targetHalfXY = 10.0 * cm;
    auto* targetSolid = new G4Box(
        "TargetSolid", targetHalfXY, targetHalfXY, 0.5 * fTargetThickness);
    auto* targetLV = new G4LogicalVolume(targetSolid, targetMaterial, "TargetLV");
    new G4PVPlacement(
        nullptr, {}, targetLV, "TargetPV", worldLV, false, 0, true);

    return worldPV;
}
