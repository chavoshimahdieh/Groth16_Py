const { ethers } = require("hardhat");
const deploy_groth16_verifier = require("./deploy_scripts/deploy_groth16_verifier_v1");
const deploy_groth16_verifier_v1 = require("./deploy_scripts/deploy_groth16_verifier_v1");

let GV_fixtureData; // Promise to store the fixture instance

async function groth16VerifierFixtureV1() {
  if (!GV_fixtureData) {
    const [owner] = await ethers.getSigners();

    
    const grInstance = await deploy_groth16_verifier_v1(owner.add);

    GV_fixtureData = {
      grInstance,
      owner,
    };
  }

  return GV_fixtureData;

  // return await arInstancePromise; // Await the promise to get the resolved value
}

module.exports = {
  groth16VerifierFixtureV1,
};
