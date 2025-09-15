// card.tsx
import * as React from "react";

type DivProps = React.HTMLAttributes<HTMLDivElement>;
type HeadingProps = React.HTMLAttributes<HTMLHeadingElement>;

export function Card({ children, className = "", ...props }: DivProps) {
  return (
    <div
    className={`bg-primary-light shadow-md rounded-lg p-4 ${className}`}

      {...props}
    >
      {children}
    </div>
  );
}

export function CardHeader({ children, className = "", ...props }: DivProps) {
  return (
    <div className={`border-b-2 pb-4 ${className}`} {...props}>
      {children}
    </div>
  );
}

export function CardTitle({ children, className = "", ...props }: HeadingProps) {
  return (
    <h2
      className={`text-xl font-semibold text-gray-800 ${className}`}
      {...props}
    >
      {children}
    </h2>
  );
}

export function CardContent({ children, className = "", ...props }: DivProps) {
  return (
    <div className={`mt-4 ${className}`} {...props}>
      {children}
    </div>
  );
}
