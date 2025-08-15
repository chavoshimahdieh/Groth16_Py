// Import necessary packages and libraries
const { expect } = require("chai");
const { loadFixture } = require("@nomicfoundation/hardhat-network-helpers");
const { groth16VerifierFixtureV2 } = require("./fixture/groth16VerifierV2.fixture.js");

describe("Groth16VerifierV2", function () {
  let grInstance, owner;

  beforeEach(async function () {
    ({ grInstance, owner } = await loadFixture(groth16VerifierFixtureV2));
  });

  describe("emulate", () => {
    it("should return true if everything is ok", async () => {
      const isVerified = await grInstance.connect(owner).emulate();
      expect(isVerified).to.be.true;
    });
  });
});
