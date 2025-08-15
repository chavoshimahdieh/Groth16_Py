const hre = require("hardhat");
const { ethers } = require("hardhat");

async function deploy_groth16_verifier() {
  const [deployer] = await ethers.getSigners();
  console.log("Deploying the Groth16 verifier contract...");

  // Deploy Groth16Verifier
  const groth16Verifier = await ethers.getContractFactory("Groth16Verifier");
  const gvInstance = await groth16Verifier.deploy(deployer.address);
  await gvInstance.deployed();

  console.log("Groth16Verifier Contract deployed to:", gvInstance.address);
  console.log("---------------------------------------------------------");
  return grInstance.address;
}

module.exports = deploy_groth16_verifier;
