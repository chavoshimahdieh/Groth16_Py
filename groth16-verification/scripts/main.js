const hre = require("hardhat");
const { ethers } = require("hardhat");

const deploy_groth16_verifier = require("./deploy/deploy_groth16_verifier");

async function main() {
  // deploy contracts
  await deploy_groth16_verifier();
 
}
main()
  .then(() => process.exit(0))
  .catch((error) => {
    console.error(error);
    process.exit(1);
  });
