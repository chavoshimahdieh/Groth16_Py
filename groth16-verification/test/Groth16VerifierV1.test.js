// Import necessary packages and libraries
const { expect } = require("chai");
const { loadFixture } = require("@nomicfoundation/hardhat-network-helpers"); 
const { groth16VerifierFixtureV1 } = require("./fixture/groth16VerifierV1.fixture.js");

describe("Groth16VerifierV1", function () {
  let grInstance, owner;

  beforeEach(async function () {
    ({ grInstance, owner } = await loadFixture(groth16VerifierFixtureV1));
  });

  describe("emulate", () => {
    it("should return true if everything is ok", async () => {
      const isVerified = await grInstance.connect(owner).emulate();
      expect(isVerified).to.be.true;
    });
  });
});
