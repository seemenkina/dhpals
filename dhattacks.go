package dhpals

import (
	"bytes"
	"crypto/rand"
	"fmt"
	"math/big"
)

func findSmallFactors(cofactor *big.Int) []*big.Int {
	var factors []*big.Int
	for i := 2; i < 1<<16; i++ {
		I := new(big.Int).SetInt64(int64(i))
		z := new(big.Int).Mod(cofactor, I)
		if z.Cmp(Big0) == 0 {
			factors = append(factors, I)
			for new(big.Int).Mod(cofactor, I).Cmp(Big0) == 0 {
				cofactor = new(big.Int).Div(cofactor, I)
			}
		}
		if I.Cmp(cofactor) == 1 {
			break
		}
	}
	return factors
}

// find an element of order order mod p
func findElementOfOrder(order, p *big.Int) (*big.Int, error) {
	h := Big1

	for h.Cmp(Big1) == 0 {
		buf, err := rand.Int(rand.Reader, p)
		if err != nil {
			return nil, err
		}
		if buf.Cmp(Big0) == 0 {
			buf = buf.Add(buf, Big1)
		}
		expBuf := new(big.Int).Sub(p, Big1)
		exp := new(big.Int).Div(expBuf, order)

		h = new(big.Int).Exp(buf, exp, p)
	}

	return h, nil
}

func bruteForcedMAC(digest []byte, h, order, p *big.Int) (*big.Int, error) {

	if new(big.Int).Exp(h, order, p).Cmp(Big1) != 0 {
		return nil, fmt.Errorf("element h should have had order order")
	}

	for i := int64(1); i <= order.Int64(); i++ {
		public := new(big.Int).Exp(h, new(big.Int).SetInt64(i), p)
		candidate := mixKey(public.Bytes())

		if bytes.Compare(candidate, digest) == 0 {
			return new(big.Int).SetInt64(i), nil
		}
	}

	return nil, fmt.Errorf(" could not find h^x such that h^x = K mod p ")
}

func runDHSmallSubgroupAttack(p, cofactor *big.Int, dh func(*big.Int) []byte) (priv *big.Int) {
	var B, R []*big.Int

	factors := findSmallFactors(cofactor)

	for _, order := range factors {
		h, err := findElementOfOrder(order, p)
		if err != nil {
			break
		}
		digest := dh(h)

		b, err := bruteForcedMAC(digest, h, order, p)
		if err != nil {
			fmt.Printf("Error: %v", err)
			break
		}
		B = append(B, b)
		R = append(R, order)
	}

	priv, _, err := crt(B, R)
	if err != nil {
		fmt.Printf("Error: %v", err)
		return nil
	}

	return priv
}

// catchKangaroo implements Pollard's kangaroo algorithm.
func catchKangaroo(p, g, y, a, b *big.Int) (m *big.Int, err error) {
	panic("not implemented")
	return
}

func runDHKangarooAttack(p, g, q, cofactor *big.Int, dh func(*big.Int) []byte, getPublicKey func() *big.Int) (priv *big.Int) {
	panic("not implemented")
	return nil
}
