{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}

-- | FR3D provides a very convenient library of explored RNA structures. We are
-- mostly interested in the "basepairs" files. In contrast to the RNAstrand
-- library or melting experiments, these data sets provide non-canonical RNA
-- pairing.

module Biobase.FR3D where

import Data.ByteString.Char8 as BS
import Data.List as L



-- | Encapsulates all the "basepairs" information.

data FR3D = FR3D
  { pdbid :: ByteString
  , chains :: [(ByteString,ByteString)]
  , basepairs :: [Basepair]
  } deriving (Show)

-- | A single basepair in a basepair system.

data Basepair = Basepair
  { interaction :: ByteString
  -- nucleotide 1
  , nucleotide1 :: Char
  , pdbnumber1 :: Int
  , chain1 :: ByteString
  , seqpos1 :: Int
  -- nucleotide 2
  , nucleotide2 :: Char
  , pdbnumber2 :: Int
  , chain2 :: ByteString
  , seqpos2 :: Int
  } deriving (Show)

-- | Linearized FR3D format.

data LinFR3D = LinFR3D
  { pdbID :: ByteString
  , sequence :: ByteString
  , pairs :: [(Int,Int,String)] -- TODO String -> CWW ?!
  } deriving (Show)

-- | The default format is a bit unwieldy; Linearization assumes that all
-- sequences are in 5'->3' order; then produces one sequence with "&"
-- separating the sequences and pairs reduced to (Int,Int,cWW).

linearizeFR3D :: FR3D -> LinFR3D
linearizeFR3D FR3D{..} = LinFR3D
  { pdbID = pdbid
  , sequence = BS.intercalate "&" $ L.map snd chains
  , pairs = L.map f basepairs
  } where
      trans = snd $ L.mapAccumL ( \acc (x,y) -> (acc + 1 + BS.length y, (x,acc))
                                ) 0 chains
      f Basepair{..} =  ( maybe (-1) (\v -> v+pdbnumber1) $ L.lookup chain1 trans
                        , maybe (-1) (\v -> v+pdbnumber2) $ L.lookup chain2 trans
                        , BS.unpack interaction
                        )
