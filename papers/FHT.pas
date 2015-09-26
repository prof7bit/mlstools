function FHT(const Data: TIntegerArray; const Poly, Seed: Cardinal): TIntegerArray; overload;

  function Log2(const Value: Cardinal): Integer;
  // Result := Ln2(Value)
  asm
         BSR   EAX,EAX
         JNZ   @@1
         DEC   EAX
  @@1:
  end;

  function Permute(const Data: TIntegerArray; Poly, Seed: Cardinal): TIntegerArray;
  // nutze inverse MLS zum Polynom "Poly" um die Eingangssamples umzusortieren für die nachfolgende FHT
  // "Poly" ist also das Polynom der MLS kodierten Daten
  var
    I,Bits,Size,DC,Value: Integer;
    Mask,Idx: Cardinal;
  begin
    Bits := Log2(Poly);
    Mask := 1 shl Bits;
    Size := Mask -1;
    Seed := Seed and Size;
    Idx  := 0;
    Poly := Poly shr 1;
    Mask := Mask shr 1;
    for I := 0 to Bits -2 do
    begin
      Idx := Idx shr 1;
      if Odd(Seed) then
      begin
        Idx  := Idx or Mask;
        Seed := (Seed shr 1) xor Poly;
      end else Seed := Seed shr 1;
    end;
    SetLength(Result, Size +1);
    DC := 0;
    for I := 0 to Size -1 do
    begin
      Idx := Idx shr 1;
      if Odd(Seed) then
      begin
        Idx  := Idx or Mask;
        Seed := (Seed shr 1) xor Poly;
      end else Seed := Seed shr 1;
      Value := Data[I];
      Dec(DC, Value);
      Result[Idx] := Value;
    end;
    Result[0] := DC;
  end;

  function InvPermute(const Data: TIntegerArray; Poly, Seed: Cardinal): TIntegerArray;
  // inverse Permutation der Ausgangsdaten der FHT
  var
    I,Bits,Size: Integer;
  begin
    Bits := Log2(Poly);
    Size := 1 shl Bits -1;
    Seed := Poly;
    Poly := 0;
    for I := 0 to Bits -1 do
    begin
      Inc(Poly, Poly + Seed and 1);
      Seed := Seed shr 1;
    end;
    SetLength(Result, Size);
    for I := 0 to Size -1 do
    begin
      if Odd(Seed) then Seed := (Seed shr 1) xor Poly
        else Seed := Seed shr 1;
      Result[I] := Data[Seed];
    end;
  end;

  function Butterfly(Data: TIntegerArray): TIntegerArray;
  // FHT selber
  var
    I,J,K1,K2,P,T: Integer;
  begin
    P  := Length(Data);
    K2 := P;
    repeat
      K1 := K2;
      K2 := K2 div 2;
      for J := 0 to K2 -1 do
      begin
        I := J;
        repeat
          T            := Data[I + K2];
          Data[I + K2] := Data[I] - T;
          Data[I]      := Data[I] + T;
          Inc(I, K1);
        until I >= P;
      end;
    until K2 = 1;
    Result := Data;
  end;

begin
  Result := InvPermute(Butterfly(Permute(Data, Poly, Seed)), Poly, Seed);
end;
