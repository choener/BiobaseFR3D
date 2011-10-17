{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards #-}

-- | FR3D provides a very convenient library of explored RNA structures. We are
-- mostly interested in the "basepairs" files. In contrast to the RNAstrand
-- library or melting experiments, these data sets provide non-canonical RNA
-- pairing.
--
-- NOTE that FR3D entries contain basepairs both in (i,j) as well as (j,i)
-- orientation (with i<j).

module Biobase.FR3D where

import Data.ByteString.Char8 as BS
import Data.List as L

import Biobase.Secondary



-- | Encapsulates all the "basepairs" information.

data FR3D = FR3D
  { pdbid :: ByteString
  , chains :: [(ByteString,ByteString)]
  , basepairs :: [Basepair]
  } deriving (Show)

-- | A single basepair in a basepair system.

data Basepair = Basepair
  { interaction :: ExtPairAnnotation
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
  , pairs :: [ExtPairIdx] -- [(Int,Int,String)] -- TODO String -> CWW ?!
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
      f Basepair{..} =  ( ( maybe (-1) (\v -> v+seqpos1) $ L.lookup chain1 trans
                          , maybe (-1) (\v -> v+seqpos2) $ L.lookup chain2 trans )
                        , interaction
                        )

class RemoveDuplicatePairs a where
  removeDuplicatePairs :: a -> a

instance RemoveDuplicatePairs FR3D where
  removeDuplicatePairs x@FR3D{..} = x{basepairs = L.filter f basepairs} where
    f Basepair{..} = (chain1,seqpos1) < (chain2,seqpos2)

instance RemoveDuplicatePairs LinFR3D where
  removeDuplicatePairs x@LinFR3D{..} = x{pairs = L.filter f pairs} where
    f ((x,y),_) = x<y


-- ** Checking data structures
--
-- Two functions to check 'FR3D' and 'LinFR3D' data structures.

-- | Checks an FR3D file for correctness. Returns either a Left on errors or
-- Right FR3D if correct.
--
-- TODO chain existence check

checkFR3D fr3d@FR3D{..}
  | L.null xs = Right fr3d
  | otherwise = Left xs
  where
    xs = [ x
         | x <- basepairs
         , let Just c1 = lookup (chain1 x) chains
         , let Just c2 = lookup (chain2 x) chains
         , nucleotide1 x /= c1 `BS.index` seqpos1 x || nucleotide2 x /= c2 `BS.index` seqpos2 x
         ]

checkLinFR3D linfr3d@LinFR3D{..}
  | L.null xs = Right linfr3d
  | otherwise = Left xs
  where
    xs = [ x
         | x <- pairs
