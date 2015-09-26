{ follow procedures create and check suitable Polynoms for above
  Linear Feedback Shift Registers.
  Such polynoms must be primitive, irreducible to produce maximal
  length LFSR's. }


type
  TFactors = array[0..6] of Cardinal;

function FactorOrder(var Factors: TFactors; Degree: Integer): Integer;
// find factors of 2^Degree-1 = possible Order of Polynom
// can be surrely more optimized, but the runtime here is not important yet
// as example:
//   For all orders from 2^0 -1 upto 2^32-1 exists only 33 possible primefactors.
//   Instead to looping trough all odd numbers as factors we could reduce the
//   loopcount if we use a lookuptable of all the 33 possible primefactors.
var
  Order,Rest,Prime,Bound,Factor,LastFactor: Cardinal;
begin
  Result := 0;
  LastFactor := 0;
  Prime := 3;
  Order := $FFFFFFFF shr (32 - Degree);
  Rest := Order;
  Bound := Round(Sqrt(Rest));
  while (Rest <> 1) and (Prime < Bound) do
    if Rest mod Prime = 0 then
    begin
      Rest := Rest div Prime;
      Factor := Order div Prime;
      if Factor <> LastFactor then
      begin
        LastFactor := Factor;
        Factors[Result] := Factor;
        Inc(Result);
      end;
      Bound := Round(Sqrt(Rest));
    end else Inc(Prime, 2);
  if Result > 0 then
  begin // only if 2^Degree-1 it self isn't prime
    Factors[Result] := Order div Rest;
    Inc(Result);
  end;
end;

function PolyIsPrimitive(Poly: Cardinal; Degree: Integer; const Factors: TFactors; FactorCount: Integer): Boolean;
// small, clean and efficient for small Degrees upto 32
// for larger Degrees we should go with real GF(2) arithmetic
type
  TPolyMatrix = array[0..31] of Cardinal;

  procedure PolyCopy(var D: TPolyMatrix; const S: TPolyMatrix);
  {$IFDEF ASM}
  asm
         XCHG  EAX,EDI
         XCHG  EDX,ESI
         MOV   ECX,32
         REP   MOVSD
         MOV   EDI,EAX
         MOV   ESI,EDX
  end;
  {$ELSE}
  begin
    Move(S, D, SizeOf(S));
  end;
  {$ENDIF}

  procedure PolyInit(var M: TPolyMatrix; Poly: Cardinal; Degree: Integer);
  var
    L: Cardinal;
    I: Integer;
  begin
    L := $80000000;
    M[Degree -1] := L;
    for I := 0 to Degree -2 do
    begin
      L    := L shr 1;
      M[I] := L;
    end;
    for I := 0 to Degree -2 do
    begin
      if Odd(Poly) then M[I] := M[I] or $80000000;
      Poly := Poly shr 1;
    end;
  end;

  procedure PolyMul(var R: TPolyMatrix; const M: TPolyMatrix; Degree: Integer);
  {$IFDEF ASM}
  // speedup >40% over all, and of course follow asm is nothing important or special tricky
  asm
         PUSH  EDI
         PUSH  ESI
         PUSH  EBX                       
         PUSH  EBP
         LEA   EDI,[EAX + ECX * 4]  // on top of R, eg @R[Degree]
         LEA   ESI,[EDX + ECX * 4]  // on top of M
         NEG   ECX                  // -Degree, we loop from -Degree upto 0
         MOV   EBP,ECX
  @@1:   MOV   EDX,[ESI + EBP * 4]  // outerloop
         XOR   EAX,EAX
         PUSH  ECX                  // save -Degree, eg. loopcounter
  @@2:   ADD   EDX,EDX              // inner loop, branchfree xor'ing is important for speed
         SBB   EBX,EBX
         AND   EBX,[EDI + ECX * 4]
         XOR   EAX,EBX
         INC   ECX
         JNZ   @@2
         POP   ECX                  // restore loopcounter
         INC   EBP
         PUSH  EAX                  // push T[] on stack
         JNZ   @@1
         SUB   EDI,4
  @@3:   POP   DWord Ptr [EDI]      // R[] = T[]
         SUB   EDI,4
         INC   ECX
         JNZ   @@3
         POP   EBP
         POP   EBX
         POP   ESI
         POP   EDI
  end;
  {$ELSE}
  var
    T: TPolyMatrix;
    I,J: Integer;
    N,D: Cardinal;
  begin
    for I := 0 to Degree -1 do
    begin
      N := M[I];
      D := 0;
      for J := 0 to Degree -1 do
      begin
        if N and $80000000 <> 0 then D := D xor R[J];
        N := N shl 1;
      end;
      T[I] := D;
    end;
    PolyCopy(R, T);
  end;
  {$ENDIF}

  procedure PolyPowMod(var R: TPolyMatrix; const M: TPolyMatrix; N: Cardinal; Degree: Integer);
  var
    L: Cardinal;
  begin
    PolyCopy(R, M);
    L := $80000000;
    while N and L = 0 do L := L shr 1;
    while L > 1 do
    begin
      L := L shr 1;
      PolyMul(R, R, Degree);
      if L and N <> 0 then PolyMul(R, M, Degree);
    end;
  end;

var
  P,M: TPolyMatrix;
  I,J: Integer;
  State: Cardinal;
begin
  Result := False;
  PolyInit(M, Poly, Degree);
  PolyCopy(P, M);
  State := P[0];
  for I := 1 to Degree do
  begin
    PolyMul(P, P, Degree);
    if P[0] = State then
      if I = Degree then
      begin
        for J := 0 to FactorCount -1 do
        begin
          PolyPowMod(P, M, Factors[J], Degree);
          if P[0] = $80000000 then Exit;
        end;
        Result := True;
        Exit;
      end else Exit;  
  end;    
end;

function FindPrimitivePolynom(Poly: Cardinal; Degree: Integer): Cardinal;
// on PIV 1.5 GHz
// 4.95 ms (asm) 6.97 ms (pascal) on degree 32 searches
// 0.41 ms (asm) 0.60 ms (pascal) on degree 16 searches
// on random input about 1/16 are primitive polynoms
// Poly / 2, means Result *2 +1 = real polynom

  function EvenParity(Poly: Cardinal): Boolean;
  // returns TRUE if count of bits set to 1 in Poly is even
  {$IFDEF ASM}
  asm
       MOVZX  EDX,AX
       SHR    EAX,16
       XOR    EAX,EDX
       XOR     AH,AL
       SETP    AL
  end;
  {$ELSE}
  var
    I: Cardinal;
  begin
    I := 0;
    while Poly <> 0 do
    begin
      Inc(I, Poly and 1);
      Poly := Poly shr 1;
    end;
    Result := not Odd(I);
  end;
  {$ENDIF}

var
  Factors: TFactors;
  Mask: Cardinal;
  FactorCount: Integer;
begin
  if Degree > 0 then
  begin
    Mask := $FFFFFFFF shr (32 - Degree);   // eg. 2^Degree -1
    Poly := Poly and Mask or 1 shl (Degree -1);
    FactorCount := FactorOrder(Factors, Degree);
    while Poly <= Mask do
      if PolyIsPrimitive(Poly, Degree, Factors, FactorCount) then
      begin
        Result := Poly;
        Exit;
      end else
      repeat
        Inc(Poly);
      until (Poly > Mask) or EvenParity(Poly);
  end;
  Result := 0;
end;

function IsPolynomPrimitive(Poly: Cardinal): Boolean;
// on PIV 1.5 GHz
// ~0.446 ms on random polys with degree 32
// ~0.052 ms on random polys with degree 16
var
  Factors: TFactors;
  FactorCount,Degree: Integer;
  P: Cardinal;
begin
  Degree := 0;
  P := Poly;
  while P <> 0 do
  begin
    P := P shr 1;
    Inc(Degree);
  end;
  if Degree > 0 then
  begin
    FactorCount := FactorOrder(Factors, Degree);
    Result := PolyIsPrimitive(Poly, Degree, Factors, FactorCount);
  end else Result := False;
end;
