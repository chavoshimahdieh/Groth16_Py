const hre = require("hardhat");
const { ethers } = require("hardhat");

async function deploy_groth16_verifier_v2(deployer) {

  // Deploy Groth16Verifier
  const groth16Verifier = await ethers.getContractFactory("Groth16VerifierV2");
  const grInstance = await groth16Verifier.deploy();
  await grInstance.waitForDeployment();

  return grInstance;
}

module.exports = deploy_groth16_verifier_v2;
