const { ethers } = require("hardhat");
const deploy_groth16_verifier_v2 = require("./deploy_scripts/deploy_groth16_verifier_v2");

let GV_fixtureData; // Promise to store the fixture instance

async function groth16VerifierFixtureV2() {
  if (!GV_fixtureData) {
    const [owner] = await ethers.getSigners();

    
    const grInstance = await deploy_groth16_verifier_v2(owner.add);

    GV_fixtureData = {
      grInstance,
      owner,
    };
  }

  return GV_fixtureData;

  // return await arInstancePromise; // Await the promise to get the resolved value
}

module.exports = {
  groth16VerifierFixtureV2,
};
